# Licensed under a 3-clause BSD style license - see LICENSE.rst
import pytest

from astropy.coordinates import EarthLocation, SkyCoord
from astropy import units as u
from astropy.time import Time

try:
    import matplotlib  # noqa
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# TODO: replace this with actual plot checks once these
# issues are resolved:
# https://github.com/astropy/astroplan/issues/65
# https://github.com/astropy/astroplan/issues/74
@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.mpl_image_compare(baseline_dir='baseline_images')
def test_image_example():
    import matplotlib.pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot([1, 2, 3])

    return fig


@pytest.mark.skipif('not HAS_MATPLOTLIB')
@pytest.mark.mpl_image_compare(baseline_dir='baseline_images')
def test_timezone():
    import datetime
    import zoneinfo

    import matplotlib.pyplot as plt

    from astroplan import Observer
    from astroplan.plots.time_dependent import plot_airmass

    betelgeuse = SkyCoord(88.79293899*u.deg, 7.407064*u.deg, frame='icrs')
    subaru = EarthLocation.from_geodetic(-155.4761111111111*u.deg,
                                         19.825555555555564*u.deg, 4139*u.m)
    observer = Observer(subaru)
    # Eastern time... because you're remote-operating Subaru from home...?
    # Use a fixed time so the image is reproducible.
    time_ET = datetime.datetime(2020, 1, 15, 20, 0,
                                tzinfo=zoneinfo.ZoneInfo('US/Eastern'))

    fig, ax = plt.subplots()
    plot_airmass(betelgeuse, observer, time_ET, ax=ax, use_local_tz=True)
    fig.tight_layout()

    return fig


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_plot_altitude():
    import matplotlib.pyplot as plt

    from astroplan import Observer
    from astroplan.target import FixedTarget
    from astroplan.plots.time_dependent import plot_altitude

    location = EarthLocation.from_geodetic(lon=0*u.deg, lat=51*u.deg, height=0*u.m)
    observer = Observer(location=location)
    target = FixedTarget(coord=SkyCoord(ra=10*u.deg, dec=45*u.deg), name='test')
    time = Time('2024-01-01 00:00:00')

    fig, ax = plt.subplots()
    result = plot_altitude(target, observer, time, ax=ax)
    assert result is ax
    plt.close(fig)


@pytest.mark.skipif('not HAS_MATPLOTLIB')
def test_plot_parallactic():
    import matplotlib.pyplot as plt

    from astroplan import Observer
    from astroplan.target import FixedTarget
    from astroplan.plots.time_dependent import plot_parallactic

    location = EarthLocation.from_geodetic(lon=0*u.deg, lat=51*u.deg, height=0*u.m)
    observer = Observer(location=location)
    target = FixedTarget(coord=SkyCoord(ra=10*u.deg, dec=45*u.deg), name='test')
    time = Time('2024-01-01 00:00:00')

    fig, ax = plt.subplots()
    result = plot_parallactic(target, observer, time, ax=ax)
    assert result is ax
    plt.close(fig)
