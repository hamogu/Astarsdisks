from astropy import units as u
from sherpa.astro import ui


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