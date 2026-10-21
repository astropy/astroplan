# Licensed under a 3-clause BSD style license - see LICENSE.rst

from unittest.mock import patch

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time

from astroplan import FixedTarget, Observer, Schedule
from astroplan.plots import plot_airmass, plot_altitude, plot_schedule_airmass

plt = pytest.importorskip("matplotlib.pyplot")


def _observer(latitude):
    return Observer(latitude=latitude * u.deg, longitude=0 * u.deg)


def _target():
    return FixedTarget(coord=SkyCoord(ra=0 * u.deg, dec=20 * u.deg), name="Target")


@pytest.mark.parametrize("plot_function", [plot_airmass, plot_altitude])
@pytest.mark.parametrize(
    "latitude, date, has_shading",
    [
        (35, "2020-06-19", True),
        (60, "2020-06-19", True),
        (75, "2020-06-19", False),
        (75, "2020-12-19", True),
    ],
)
def test_night_shading_latitudes(plot_function, latitude, date, has_shading):
    observer = _observer(latitude)
    times = Time(date) + np.linspace(0, 3, 13) * u.day
    fig, ax = plt.subplots()

    plot_function(_target(), observer, times, ax=ax, brightness_shading=True)

    assert bool(ax.patches) == has_shading
    assert all(shading_patch.get_alpha() > 0 for shading_patch in ax.patches)

    if latitude == 60:
        # Astronomical twilight does not occur in summer.
        assert all(shading_patch.get_alpha() < 0.4 for shading_patch in ax.patches)

    plt.close(fig)


def test_night_shading_does_not_cover_daytime():
    observer = _observer(49)
    times = Time("2020-06-19") + np.linspace(0, 4, 13) * u.day
    noon = Time("2020-06-20 12:00:00").plot_date
    fig, ax = plt.subplots()

    plot_airmass(_target(), observer, times, ax=ax, brightness_shading=True)

    assert ax.patches

    for shading_patch in ax.patches:
        x = (
            shading_patch.get_path()
            .transformed(shading_patch.get_patch_transform())
            .vertices[:, 0]
        )
        assert not (x.min() <= noon <= x.max())

    plt.close(fig)


def test_night_shading_vectorized():
    observer = _observer(35)
    times = Time("2020-06-19") + np.linspace(0, 5, 13) * u.day
    fig, ax = plt.subplots()

    with (
        patch.object(observer, "sun_rise_time", wraps=observer.sun_rise_time) as rise,
        patch.object(observer, "sun_set_time", wraps=observer.sun_set_time) as set_,
    ):
        plot_airmass(_target(), observer, times, ax=ax, brightness_shading=True)

    assert rise.call_count == 4
    assert set_.call_count == 4

    plt.close(fig)


def test_empty_schedule_night_shading():
    schedule = Schedule(Time("2020-12-19"), Time("2020-12-22"))
    schedule.observer = _observer(35)
    fig, ax = plt.subplots()
    plt.sca(ax)

    plot_schedule_airmass(schedule, show_night=True)

    assert ax.patches
    np.testing.assert_allclose(
        ax.get_xlim(), [schedule.start_time.plot_date, schedule.end_time.plot_date]
    )

    plt.close(fig)


@pytest.mark.mpl_image_compare(baseline_dir="baseline_images")
def test_multi_night_shading_image():
    observer = _observer(35)
    times = Time("2020-06-19") + np.linspace(0, 3, 13) * u.day
    fig, ax = plt.subplots()

    plot_airmass(_target(), observer, times, ax=ax, brightness_shading=True)

    return fig
