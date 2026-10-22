# Licensed under a 3-clause BSD style license - see LICENSE.rst

from datetime import datetime, timedelta, timezone

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time

from astroplan import FixedTarget, Observer, Schedule
from astroplan.plots import plot_airmass, plot_altitude, plot_schedule_airmass
from astroplan.plots.time_dependent import _shade_nights

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
    times = Time(date) + np.linspace(0, 3, 5) * u.day
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


@pytest.mark.parametrize("date", [
    "2020-03-24",
    "2020-04-14",
    "2020-05-14",
])
def test_night_shading_boundaries(date):
    observer = Observer(latitude=32 * u.deg, longitude=-111 * u.deg)
    start = Time(date)
    end = start + 1 * u.day
    fig, ax = plt.subplots()

    _shade_nights(ax, observer, start, end)

    # Find expected crossings with a one-minute altitude scan.
    scan_times = start + np.arange(24 * 60 + 1) * u.min
    sun_altitudes = observer.sun_altaz(scan_times).alt.to_value(u.deg)
    bands = np.searchsorted([-18, -12, -6, 0], sun_altitudes)
    expected_edges = scan_times.plot_date[1:][np.diff(bands) != 0]

    # Extract the boundaries of the plotted shading.
    edges = []
    for shading_patch in ax.patches:
        x = shading_patch.get_path().transformed(
            shading_patch.get_patch_transform()).vertices[:, 0]
        edges.extend([x.min(), x.max()])

    edges = np.unique(edges)
    edges = edges[(edges > start.plot_date + 1e-7) &
                  (edges < end.plot_date - 1e-7)]

    np.testing.assert_allclose(
        edges, expected_edges, rtol=0, atol=2 / (24 * 60)
    )

    plt.close(fig)


def test_empty_schedule_night_shading():
    schedule = Schedule(Time("2020-12-19"), Time("2020-12-22"))
    schedule.observer = _observer(35)
    fig, ax = plt.subplots()
    plt.sca(ax)

    plot_schedule_airmass(schedule, show_night=True)

    assert ax.patches
    assert ax.get_ylabel() == "Airmass"
    assert ax.get_xlabel() == "Time from 2020-12-19 [UTC]"
    np.testing.assert_allclose(ax.get_ylim(), [3.0, 1.0])
    np.testing.assert_allclose(
        ax.get_xlim(),
        [schedule.start_time.plot_date, schedule.end_time.plot_date],
        rtol=0, atol=1e-8
    )

    plt.close(fig)


def test_night_shading_local_timezone():
    observer = _observer(35)
    time = datetime(2020, 6, 19, 20, tzinfo=timezone(timedelta(hours=2)))
    fig, (ax_utc, ax_local) = plt.subplots(2, 1)

    plot_airmass(_target(), observer, time, ax=ax_utc,
                 brightness_shading=True)
    plot_airmass(_target(), observer, time, ax=ax_local,
                 brightness_shading=True, use_local_tz=True)

    offset = 2 / 24

    np.testing.assert_allclose(
        ax_local.lines[0].get_xdata(),
        ax_utc.lines[0].get_xdata() + offset,
        rtol=0, atol=1e-8
    )

    assert ax_utc.patches
    assert len(ax_local.patches) == len(ax_utc.patches)

    for utc_patch, local_patch in zip(ax_utc.patches, ax_local.patches):
        utc_x = utc_patch.get_path().transformed(
            utc_patch.get_patch_transform()).vertices[:, 0]
        local_x = local_patch.get_path().transformed(
            local_patch.get_patch_transform()).vertices[:, 0]

        np.testing.assert_allclose(
            local_x, utc_x + offset, rtol=0, atol=1e-8
        )
        assert local_patch.get_alpha() == utc_patch.get_alpha()

    plt.close(fig)


@pytest.mark.parametrize("plot_function", [plot_airmass, plot_altitude])
def test_night_shading_descending_times(plot_function):
    observer = _observer(35)
    ascending = Time("2020-06-19") + np.linspace(0, 3, 5) * u.day
    descending = ascending[::-1]
    fig, (ax_ascending, ax_descending) = plt.subplots(2, 1)

    plot_function(_target(), observer, ascending, ax=ax_ascending,
                  brightness_shading=True)
    plot_function(_target(), observer, descending, ax=ax_descending,
                  brightness_shading=True)

    assert ax_ascending.patches
    assert len(ax_ascending.patches) == len(ax_descending.patches)

    for ascending_patch, descending_patch in zip(
        ax_ascending.patches, ax_descending.patches
    ):
        ascending_x = ascending_patch.get_path().transformed(
            ascending_patch.get_patch_transform()).vertices[:, 0]
        descending_x = descending_patch.get_path().transformed(
            descending_patch.get_patch_transform()).vertices[:, 0]

        np.testing.assert_allclose(
            [ascending_x.min(), ascending_x.max()],
            [descending_x.min(), descending_x.max()],
            rtol=0, atol=1e-8
        )
        assert ascending_patch.get_alpha() == descending_patch.get_alpha()

    plt.close(fig)


@pytest.mark.mpl_image_compare(baseline_dir="baseline_images")
def test_multi_night_shading_image():
    observer = _observer(35)
    times = Time("2020-06-19") + np.linspace(0, 3, 289) * u.day
    fig, ax = plt.subplots()

    plot_airmass(_target(), observer, times, ax=ax, brightness_shading=True)

    fig.tight_layout()
    return fig
