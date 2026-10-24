# Licensed under a 3-clause BSD style license - see LICENSE.rst

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.time import Time

from astroplan import (AirmassConstraint, AltitudeConstraint, AtNightConstraint, Constraint,
                       FixedTarget, MoonSeparationConstraint, Observer)
from astroplan.plots import plot_airmass, plot_altitude
from astroplan.plots.time_dependent import _MET_LINEWIDTH_FACTOR

plt = pytest.importorskip("matplotlib.pyplot")

PLOT_FUNCTIONS = [plot_airmass, plot_altitude]


class _FailAt(Constraint):
    """Constraint that fails only at the given times."""
    def __init__(self, bad_times, boolean=True):
        self.bad_jd = bad_times.jd
        self.boolean = boolean

    def compute_constraint(self, times, observer, targets):
        bad = np.isclose(times.jd[..., np.newaxis], self.bad_jd, rtol=0, atol=1e-8)
        met = ~np.any(bad, axis=-1)
        return met if self.boolean else met.astype(float)


@pytest.fixture
def observer():
    return Observer(latitude=0 * u.deg, longitude=0 * u.deg)


def _target(dec, name, ra=0):
    return FixedTarget(coord=SkyCoord(ra=ra * u.deg, dec=dec * u.deg), name=name)


@pytest.fixture
def times(observer):
    """Thirteen samples spanning +/-3 h around the transit of an equatorial target."""
    transit = observer.target_meridian_transit_time(Time('2020-06-19'), _target(0, 'T'),
                                                    which='next')
    return transit + np.linspace(-3, 3, 13) * u.hour


def _curves(ax):
    """Split the lines into (curves where constraints are met, dashed curves)."""
    dashed = [line for line in ax.lines if line.get_label() == '_nolegend_']
    met = [line for line in ax.lines if line.get_label() != '_nolegend_']
    return met, dashed


def _finite_x(line):
    """Times at which the line is actually drawn (y is not NaN)."""
    x, y = np.asarray(line.get_xdata()), np.asarray(line.get_ydata(), dtype=float)
    return x[np.isfinite(y)]


def _extend(mask):
    """Extend ``mask`` by one sample on each side."""
    return np.convolve(mask, [1, 1, 1], mode='same') > 0


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
def test_no_constraints_no_dashed_curve(plot_function, observer, times):
    fig, ax = plt.subplots()
    plot_function(_target(0, 'T'), observer, times, ax=ax)
    met, dashed = _curves(ax)
    assert len(met) == 1 and not dashed
    np.testing.assert_allclose(_finite_x(met[0]), times.plot_date)
    # Without constraints the line width is untouched
    assert met[0].get_linewidth() == 1.5
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
def test_curve_is_dashed_where_constraints_fail(plot_function, observer, times):
    target = _target(0, 'T')
    fig, ax = plt.subplots()
    plot_function(target, observer, times, ax=ax,
                  constraints=[AltitudeConstraint(min=70 * u.deg)])

    failing = (observer.altaz(times, target).alt < 70 * u.deg)
    assert failing.any() and not failing.all()

    (met,), (dashed,) = _curves(ax)
    np.testing.assert_allclose(_finite_x(met), times.plot_date[~failing])
    np.testing.assert_allclose(_finite_x(dashed), times.plot_date[_extend(failing)])
    assert dashed.get_linestyle() == '--'
    assert met.get_linestyle() == '-'
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
@pytest.mark.parametrize("style_kwargs, linewidth", [
    (None, 1.5),
    ({'lw': 3}, 3),
    ({'linewidth': 0.5}, 0.5),
    ({'color': 'tab:red', 'alpha': 0.8}, 1.5),
    ({'c': 'g'}, 1.5),
    ({'linestyle': ':'}, 1.5),
    ({'ls': ':'}, 1.5),
])
def test_dashed_curve_style(plot_function, style_kwargs, linewidth, observer, times):
    fig, ax = plt.subplots()
    plot_function(_target(0, 'T'), observer, times, ax=ax, style_kwargs=style_kwargs,
                  constraints=[AltitudeConstraint(min=70 * u.deg)])

    (met,), (dashed,) = _curves(ax)
    assert met.get_label() == 'T'
    assert met.get_color() == dashed.get_color()
    assert dashed.get_alpha() == met.get_alpha()
    # Dashed stretches keep the requested width, the others are thicker
    assert dashed.get_linewidth() == pytest.approx(linewidth)
    assert met.get_linewidth() == pytest.approx(_MET_LINEWIDTH_FACTOR * linewidth)
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
def test_targets_keep_distinct_colors(plot_function, observer, times):
    fig, ax = plt.subplots()
    plot_function([_target(0, 'A'), _target(0, 'B')], observer, times, ax=ax,
                  constraints=[AltitudeConstraint(min=70 * u.deg)])

    met, dashed = _curves(ax)
    assert met[0].get_color() != met[1].get_color()
    assert [line.get_color() for line in dashed] == [line.get_color() for line in met]
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
@pytest.mark.parametrize("boolean", [True, False])
def test_single_sample_violation_is_visible(plot_function, boolean, observer, times):
    fig, ax = plt.subplots()
    plot_function(_target(0, 'T'), observer, times, ax=ax,
                  constraints=_FailAt(times[6], boolean=boolean))

    (met,), (dashed,) = _curves(ax)
    # The failing sample is dashed together with one neighbour on each side
    np.testing.assert_allclose(_finite_x(dashed), times.plot_date[5:8])
    assert times.plot_date[6] not in _finite_x(met)
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
@pytest.mark.parametrize("constraints", [
    AtNightConstraint.twilight_astronomical(),
    MoonSeparationConstraint(min=10 * u.deg),
    AirmassConstraint(max=1.5),
    [AltitudeConstraint(min=10 * u.deg), AtNightConstraint()],
    np.array([AltitudeConstraint(min=10 * u.deg), AtNightConstraint()]),
])
def test_accepts_any_constraint(plot_function, constraints, observer, times):
    fig, ax = plt.subplots()
    plot_function(_target(0, 'T'), observer, times, ax=ax, constraints=constraints)
    met, dashed = _curves(ax)
    assert len(met) == 1 and len(dashed) == 1
    plt.close(fig)


@pytest.mark.parametrize("plot_function", PLOT_FUNCTIONS)
def test_each_target_uses_its_own_constraint_row(plot_function, observer, times):
    # The second target transits 3 h later, so it is high at different times
    # and each target must get its own row of the constraint mask.
    targets = [_target(0, 'A'), _target(0, 'B', ra=45)]
    fig, ax = plt.subplots()
    plot_function(targets, observer, times, ax=ax,
                  constraints=[AltitudeConstraint(min=60 * u.deg)])

    met, dashed = _curves(ax)
    assert len(met) == len(dashed) == 2
    for target, met_line, dashed_line in zip(targets, met, dashed):
        alt = observer.altaz(times, target).alt
        failing = alt < 60 * u.deg
        assert failing.any() and not failing.all()
        above_horizon = alt > 0 * u.deg
        np.testing.assert_allclose(_finite_x(met_line),
                                   times.plot_date[~failing & above_horizon])
        np.testing.assert_allclose(_finite_x(dashed_line),
                                   times.plot_date[_extend(failing) & above_horizon])
    plt.close(fig)
