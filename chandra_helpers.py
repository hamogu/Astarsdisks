import glob
import os
from contextlib import chdir
from tempfile import NamedTemporaryFile

from astropy import table
from astropy.io import fits
from astropy import units as u
import numpy as np

from ciao_contrib.cda import search
from coords.chandra import cel_to_chandra
from ciao_contrib.region.check_fov import FOVFiles
from ciao_contrib import runtool as rt
from ciao_contrib.psf_contrib import psfSize
from ciao_contrib.cda.data import download_chandra_obsids
from sherpa. astro import ui


def search_chandra_archive(coords, cache_file):
    """Search the Chandra archive for observations of the given targets.

    Parameters
    ----------
    coords : astropy.table.Table
        Table of coordinates for the targets.
    cache_file : str
        If this file exists, it is read instead of querying the Chandra archive.
        If it does not exist, is is created after querying the archive.
    """
    try:
        # I already have a local copy of the archive search.
        all_found = table.Table.read(cache_file)
    except FileNotFoundError:
        all_found = []

        for row in coords:
            coo = row['coord']
            found = search.search_chandra_archive(coo.ra.deg, coo.dec.deg)
            found = table.Table(found)
            if len(found) > 0:
                found.keep_columns(['ObsId', 'obs_date', 'RA', 'Dec', 'Instrument', 'Exposure', 'Grating'])
                found['object'] = row['object']
                all_found.append(found)
        # Convert into astropy Table for easier indexing and remove non-unique lines
        all_found = table.Table(np.unique(table.vstack(all_found)))
        all_found.sort(['object', 'ObsId'])
        all_found.write(cache_file)

    return all_found


def download_repro(datapath, obsids):
    """Download and reprocess observations.

    Parameters
    ----------
    datapath : str
        Path to the data directory, e.g. '/data/'
    obsids : list of str or int
        List of ObsIDs to download and reprocess.
    """
    needs_download = []
    for obsid in obsids:
        if not os.path.exists(os.path.join(datapath, str(obsid))):
            needs_download.append(obsid)
    if len(needs_download) > 0:
        with chdir(datapath):
            download_chandra_obsids(needs_download)

    needs_repro = []
    for obsid in np.unique(obsids):
        if not os.path.exists(os.path.join(datapath, str(obsid), 'repro')):
            needs_repro.append(str(obsid))
    if len(needs_repro) > 0:
        with chdir(datapath):
            rt.chandra_repro(indir=','.join(needs_repro), outdir="")


def chan_coords(datapath, coord, returnvals=['theta']):
    """Get the off-axis angle for given RA/DEC coordinate.

    Parameters
    ----------
    datapath : str
        Path to the data directory, e.g. '/data/obsid'
    coord : astropy.coordinates.SkyCoord
        Coordinate of the target.

    Returns
    -------
    theta : float
        Off-axis angle in arcmin.
    """
    evt2 = glob.glob(os.path.join(datapath, '*','*evt2*'))[0]
    header = fits.getheader(evt2, 'EVENTS')
    dhead = dict(header)
    # For full precision, the header would need DY_AVG, DZ_AVG, and DTH_AVG keywords.
    # However, that's missing for data that's not been recently reprocessed.
    # We don't need full precision here, at this point, we just want to know if this data
    # is usable for us at all. We can re-process all that we want to use later.
    # However, those keywords still need to be present in the header for a call to
    # cel_to_chandra, see https://github.com/cxcsds/ciao-contrib/pull/621
    for k in ["DY_AVG", "DZ_AVG", "DTH_AVG" ]:
        if k not in dhead:
            dhead[k] = 0.
    chan_coos = cel_to_chandra(dhead, coord.ra.deg, coord.dec.deg)

    return [chan_coos[v][0] for v in returnvals]


def remove_outside_of_FOV(dataroot, chandra_obs):
    """Remove chandra observations where the source is outside of the FOV from list.

    Parameters
    ----------
    dataroot : str
        Path to the data directory, e.g. '/data/' (the OBDID is added by
        this function)
    chandra_obs : astropy.table.Table
        Table of observations with ObsID and coord columns.
    """
    chandra_obs['in_FOV'] = False
    for row in chandra_obs:
        fov = glob.glob(os.path.join(dataroot, str(row['ObsId']), 'primary', '*fov1*'))[0]
        myobs = FOVFiles(fov)
        ii = myobs.inside(row['coord'].ra.deg, row['coord'].dec.deg)
        if len(ii) > 0:
            row['in_FOV'] = True

    chandra_obs = chandra_obs[chandra_obs['in_FOV']]
    chandra_obs.remove_column('in_FOV')


def set_ardlib_bpix(datapath):
    """Set badpix files in ardlib based on the header info from an evt2file

    Parameters
    ----------
    datapath : str
        Path to the data directory, e.g. '/data/obsid'
    """
    rt.ardlib.punlearn()
    fileevt2 = glob.glob(os.path.join(datapath, 'repro', '*evt2*'))[0]
    bpix = glob.glob(os.path.join(datapath, 'repro', '*bpix*'))[0]

    if fits.getval(fileevt2, 'INSTRUME') == 'HRC':
        rt.ardlib.AXAF_HRC_I_BADPIX_FILE = bpix
    else: # ACIS
        # set badpixel files for all CCDs that are used in the observation
        for i in fits.getval(fileevt2, 'DETNAM', 'EVENTS')[5:]:
            setattr(rt.ardlib, f'AXAF_ACIS{i}_BADPIX_FILE', f"{bpix}[BADPIX{i}]")
    rt.ardlib.write_params()


def psf_of_coord(datapath, coord, psffrac=.9, energy=2.3):
    """Get the PSF size for given RA/DEC coordinate.

    Parameters
    ----------
    datapath : str
        Path to the data directory, e.g. '/data/obsid'
    coord : astropy.coordinates.SkyCoord
        Coordinate of the target.
    psffrac : float
        Fraction of PSF enclosed by the radius returned.
    energy : float
        Energy in keV.

    Returns
    -------
    psf_size : float
        PSF size in arcsec.
    """
    evt2 = glob.glob(os.path.join(datapath, 'repro', '*evt2*'))[0]
    header = fits.getheader(evt2, 'EVENTS')
    dhead = dict(header)
    chan_coos = cel_to_chandra(dhead, coord.ra.deg, coord.dec.deg)
    return psfSize(energy, chan_coos['theta'][0], chan_coos['phi'][0], psffrac)


def cutout_image(datapath, coord, outfile, acis_energy='300:4000'):
    """Cut out a small region from an evt list and bin into an image.

    Parameters
    ----------
    datapath : str
        Path to the data directory, e.g. '/data/obsid'
    coord : astropy.coordinates.SkyCoord
        Coordinate of the target.
    outfile : str
        Path and filename to the output image (should end on ".fits")
    acis_energy : str
        Energy filter for ACIS, default is '300:4000' (in eV)
    """
    evt2 = glob.glob(os.path.join(datapath, 'repro', '*evt2*'))[0]
    if fits.getval(evt2, 'INSTRUME') == 'HRC':
        binning = "[bin x=::3,y=::3]"
    else: # ACIS
        binning = f"[bin x=::1,y=::1][energy={acis_energy}]"  # add energy filter for ACIS

    rt.dmcopy(infile=f"{evt2}[(x,y)=box({coord.ra.deg}d,{coord.dec.deg}d,.5',.5')]{binning}",
              outfile=outfile, option='image', clobber=True)


def dmextract_tab_and_lc(dataroot, target, bkg_by_hand,
                          acis_energy='300:4000', timebin=2000):
    dataroot_obsid = os.path.join(dataroot, str(target['ObsId']))
    obsid = target['ObsId']
    evt2 = glob.glob(os.path.join(dataroot_obsid, 'repro','*evt2*'))[0]
    set_ardlib_bpix(dataroot_obsid)
    coord = target['coord']
    psf = target['psf_radius']
    if fits.getval(evt2, 'INSTRUME') == 'HRC':
        energy_filter = ""
        dtffile = glob.glob(os.path.join(dataroot_obsid, 'primary','*dtf*'))[0]
        ccdid = ""
    else: # ACIS
        energy_filter = f"[energy={acis_energy}]"
        dtffile = None
        ccdid = chan_coords(dataroot_obsid, coord, returnvals=['chip_id'])[0]
        ccdid = f"ccd_id={ccdid},"
    rt.dmextract.punlearn()
    outfile = f'{dataroot}/{obsid}/dmextract_{target["object"].replace(" ", "_")}.fits'
    rt.dmextract(infile=f"{evt2}{energy_filter}[bin sky=circle({coord.ra.deg}d,{coord.dec.deg}d,{psf}'')]",
                 bkg=f'{evt2}{energy_filter}[bin sky={bkg_by_hand[(obsid, target["object"])]}]',
                 outfile=outfile, clobber=True)
    outlc = f'{dataroot}/{obsid}/lc_{target["object"].replace(" ", "_")}.fits'
    rt.dmextract(infile=f"{evt2}{energy_filter}[{ccdid}sky=circle({coord.ra.deg}d,{coord.dec.deg}d,{psf}'')][bin time=::{timebin}]",
                 bkg=f'{evt2}{energy_filter}[{ccdid}sky={bkg_by_hand[(obsid, target["object"])]}]',
                 outfile=outlc, opt='ltc1', exp=dtffile, clobber=True)
    return table.Table.read(outfile, hdu=1)


def specextract(dataroot, target, bkg_by_hand):
    """Run specextract on a target.

    Parameters
    ----------
    dataroot : str
        Path to the data directory, e.g. '/data/'
        (without te "obsid" part, that is added by this function)
    target : astropy.table.Row
        Row of the table of targets.
    bkg_by_hand : dict
        Dictionary of background regions by (obsid, object) tuple.
    """
    obsid = target['ObsId']
    fileevt2 = glob.glob(os.path.join(dataroot, str(obsid), 'repro','*evt2*'))[0]
    print(f'specextract for ObsID {obsid}')
    coord = target['coord']
    psf = target['psf_radius']

    set_ardlib_bpix(os.path.join(dataroot, str(obsid)))
    rt.specextract.punlearn()
    rt.specextract(infile=f"{fileevt2}[sky=circle({coord.ra.deg}d,{coord.dec.deg}d,{psf}'')]",
                   # refcoord needed if there are few (or zero) counts
                   refcoord=f"{coord.ra.deg} {coord.dec.deg}",
               outroot=f'{dataroot}/{obsid}/spec_{target["object"].replace(" ", "_")}_',
               bkgfile=f"{fileevt2}[sky={bkg_by_hand[(obsid, target['object'])]}]",
               grouptype="NONE", binspec="NONE", weight=False, correctpsf=True, clobber=True
        )


def _read_aplimits_out(filename):
    with open(filename, 'r') as f:
        line1 = f.readline()
        line2 = f.readline()
    return float(line1.split(',')[3]), float(line2.split(',')[3])

def run_ap_limits(obs):
    obs['Detection threshold'] = 0
    obs['Upper limit'] = 0.

    for row in obs:
        with NamedTemporaryFile() as tmp:
            rt.aplimits (prob_false_detection=.10, prob_missed_detection=.5,
                         bkg_rate=None, outfile=tmp.name, verbose=1, clobber=True,
                         m=row['BG_COUNTS'],
                         A_s=row['AREA'], A_b=row['BG_AREA'],
                         T_s=row['Exposure'] * 1e3, T_b=row['Exposure'] * 1e3,
                         )
            threshold, upper_limit = _read_aplimits_out(tmp.name)
        row['Detection threshold'] = threshold
        row['Upper limit'] = upper_limit
    obs['Upper limit'].unit = u.count/u.s


def count_rate_to_flux(obsid, model, count_rate, lo=0.4, hi=7.):
    """Run Sherpa fakeit for a given model and obsid and get energy flux

    Parameters
    ----------
    obsid : int or str
        Sherpa dataset identifier. The ARF and RMF is taken from this dataset, so
        it must be loaded into Sherpa before calling this function.
        This dataset is **replaced** with face data.
    model : `~sherpa.models.model.ArithmeticModel`
        Sherpa model (with parameters set) to use for fakeit.
    count_rate : float
        Count rate to be converted into a flux
    lo : float
        Lower energy limit for the flux calculation
    hi : float
        Upper energy limit for the flux calculation

    Returns
    -------
    flux : `~astropy.units.Quantity`
        Energy flux in erg/cm^2/s
    """
    ui.set_source(obsid, model)
    exposure = 1e4
    ui.fake_pha(obsid, arf=ui.get_arf(obsid), rmf=ui.get_rmf(obsid), exposure=exposure)
    # check if we have enough counts so that we are not dominated by Poisson noise
    if ui.get_data(obsid).counts.sum() < 1000:
        exposure = exposure * 10000 / ui.get_data(obsid).counts.sum()
        ui.fake_pha(obsid, arf=ui.get_arf(obsid), rmf=ui.get_rmf(obsid), exposure=exposure)
    scale = ui.get_data(obsid).counts.sum() / exposure / count_rate
    return ui.calc_energy_flux(lo=lo, hi=hi, id=obsid) / scale * u.erg/u.cm**2/u.s