from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time

from astroquery.simbad import Simbad


def coords_from_simbad(targets, cache_file):
    """Get coordinates, parallax, and proper motion from Simbad.

    Parameters
    ----------
    targets : list of str
        List of targets to query with SIMBAD identifiers or coordinates.
    cache_file : str
        Path to cache file to store results. If the file exists, it will be
        read instead of querying SIMBAD. If it does not exist, it will be
        created after querying SIMBAD.
    """
    try:
        # I already have a local copy of the Simbad query.
        coords = Table.read(cache_file)
    except FileNotFoundError:
        customSimbad = Simbad()

        customSimbad.add_votable_fields('typed_id')
        customSimbad.add_votable_fields('plx')
        customSimbad.add_votable_fields('propermotions')
        customSimbad.add_votable_fields('ra(d;A;ICRS;J2000;2000)', 'dec(d;D;ICRS;J2000;2000)')
        coords = customSimbad.query_objects(targets)
        coords.rename_column('TYPED_ID', 'object')
        coords['object'] = [str(o) for o in coords['object']]
        coords.write(cache_file)

    # Then add coordinates as astropy objects
    coords['coord'] = SkyCoord(ra=coords['RA_d_A_ICRS_J2000_2000'],
                           dec=coords['DEC_d_D_ICRS_J2000_2000'],
                           pm_ra_cosdec=coords['PMRA'],
                           pm_dec=coords['PMDEC'],
                           distance=coords['PLX_VALUE'].to(u.pc, equivalencies=u.parallax()),
                           obstime=Time('J2000', format='jyear_str', scale='tcb'),
                           frame="icrs")
    return coords


def coords_for_obs(obs, coords):
    """Set coordinates ofr each target, considering proper motions.

    Parameters
    ----------
    obs : astropy.table.Table
        Table of observations.
    coords : astropy.table.Table
        Table of coordinates for the targets.

    """
    coos = []
    for row in obs:
        coo = coords['coord'][(coords['object'] == row['object']).nonzero()[0][0]]
        coo = coo.apply_space_motion(new_obstime=Time(row['obs_date']))
        coos.append(coo)
    obs['coord'] = coos