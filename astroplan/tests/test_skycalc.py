from astropy.tests.helper import raises
from astroplan import skymodel
from astroplan.exceptions import SkyCalcError


@raises(ValueError)
def test_null_observatory():
    """Tests that skymodel raises an error if an observatory
    which SkyCalc doesn't support is given."""
    return skymodel(observatory='apo')


@raises(SkyCalcError)
def test_null_parameter():
    """Tests that skymodel raises an error if given
    a parameter which SkyCalc doesn't support."""
    return skymodel(airmass='test')


def test_skymodel_units():
    """Tests that skymodel returns the atmospheric
    transmission function as a `~astropy.units.quantity.Quantity`
    """
    sm = skymodel()
    assert hasattr(sm[0], "unit") and hasattr(sm[1], "unit")
