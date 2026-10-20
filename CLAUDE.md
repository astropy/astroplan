# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

astroplan is an Astropy-affiliated observation-planning package (rise/set times, observability constraints, scheduling, plotting). Python >= 3.12; core deps are numpy, astropy (>= 6.0), and tzdata (timezones use stdlib `zoneinfo`; pytz was dropped, but user-supplied pytz tzinfos must keep working). matplotlib and astroquery are optional (`[all]` extra).

## Commands

```bash
pip install -e ".[dev]"               # dev install: [all] + [test] + tox/tox-uv (version comes from setuptools_scm -> astroplan/_version.py)

pytest astroplan docs                 # tests + RST doctests (pyproject sets --doctest-rst)
pytest astroplan/tests/test_observer.py::test_sunrise_sunset_equator   # single test
pytest --remote-data=any              # also run @pytest.mark.remote_data tests (network)

flake8 astroplan --count              # lint: max line length 100, __init__.py excluded (.flake8)

tox -e py314-test-alldeps             # CI-style run (runs from .tmp/, installs package, MPLBACKEND=agg)
tox -e py314-test-image               # image comparison tests only (--mpl, pinned matplotlib)
tox -e build_docs                     # sphinx-build -W (warnings are errors); needs graphviz `dot`
tox -e codestyle
```

tox uses `tox-uv` (auto-provisioned via `requires` in `tox.ini`), so environments are built with uv and have no `pip`. The `oldestdeps` factor installs every dependency at the lowest version allowed (`uv_resolution = lowest`): direct lower bounds come from `pyproject.toml`, transitive ones from `oldestdeps-constraints.txt` (via `UV_CONSTRAINT`). If oldestdeps breaks, raise the relevant lower bound rather than adding exact pins. Derive transitive bounds from a cross-platform resolution (`uv pip compile --universal --resolution lowest ...`): some dependencies (e.g. keyring's `secretstorage`/`cryptography`) only install on Linux or Windows, and a single-platform resolution misses them.

Test configuration notes (from `pyproject.toml`):
- `filterwarnings = error`: **any new warning fails the test suite**. Fix the warning or add a narrowly scoped ignore.
- Docs `.rst` files under `docs/` are doctested, so code examples in tutorials must actually run. Exceptions: files marked `.. doctest-skip-all` (e.g. `docs/tutorials/periodic.rst`, whose printed outputs can go stale) and `.. plot::` directives, which only run in `tox -e build_docs`.
- `xfail_strict = true`.
- Plot tests in `astroplan/plots/tests/` use `@pytest.mark.mpl_image_compare` with baselines in `baseline_images/`; the images are only compared when pytest is run with `--mpl` (CI does this in the `py314-test-image` tox env, which pins matplotlib to the version the baselines were made with; failures upload a `results/` diff summary). Regenerate with `pytest <test> --mpl-generate-path=astroplan/plots/tests/baseline_images` using that pinned matplotlib, visually check the result before committing, and bump the pin in `tox.ini` whenever baselines are regenerated with a newer matplotlib.
- Tests needing the network (e.g. `FixedTarget.from_name`, `Observer.at_site` for non-builtin sites) must be marked `@pytest.mark.remote_data`.

## Architecture

All public names are re-exported flat from `astroplan/__init__.py` via `from .module import *`, so every module must maintain an accurate `__all__`. Importing `astroplan` calls `download_IERS_A()` as a side effect (falls back with `OldEarthOrientationDataWarning` if offline; see `utils.py` and `docs/faq/iers.rst`).

- **`observer.py` — `Observer`**: the central object. Wraps an `EarthLocation` plus timezone/weather (pressure, temperature, humidity). Provides `altaz`, `parallactic_angle`, rise/set/transit/twilight times, moon/sun helpers, `is_night`, etc. Rise/set/transit are computed on a grid over 24 h (`_generate_24hr_grid`, `n_grid_points`) and refined by interpolation of horizon crossings (`_horiz_cross`, `_two_point_interp`), with `which='next'|'previous'|'nearest'`. Many methods take `grid_times_targets=True` to broadcast times × targets into a 2D result.
- **`target.py`**: `FixedTarget` (wraps a `SkyCoord`, `from_name` uses Sesame/network) and `get_skycoord`, which normalizes any mix of `FixedTarget`/`SkyCoord`/lists into a single `SkyCoord` — most APIs accept any of these via this function.
- **`constraints.py`**: `Constraint` subclasses implement `compute_constraint(times, observer, targets)` returning either booleans or a 0–1 score (`boolean_constraint=False`, using `min_best_rescale`/`max_best_rescale`). `Constraint.__call__` handles time grids and target normalization. Expensive alt/az and moon computations are cached **on the Observer instance** (`observer._altaz_cache`, `_moon_cache`, `_meridian_transit_cache`) keyed by `_make_cache_key(times, targets)`. Top-level helpers: `is_observable`, `is_always_observable`, `is_event_observable`, `months_observable`, `observability_table`.
- **`scheduling.py`**: `ObservingBlock` (target + duration + priority + per-block constraints) is scored by `Scorer` (combines constraint scores into a time×block array). `Scheduler` subclasses (`SequentialScheduler`, `PriorityScheduler`) implement `_make_schedule` and fill a `Schedule` composed of `Slot`s; `Transitioner` inserts `TransitionBlock`s for slews and instrument reconfigurations.
- **`periodic.py`**: `PeriodicEvent` / `EclipsingSystem` for transits/eclipses; used by `PhaseConstraint`, `PrimaryEclipseConstraint`, `SecondaryEclipseConstraint`, and `is_event_observable`.
- **`moon.py`**: moon phase angle / illumination.
- **`plots/`**: matplotlib-based `plot_airmass`, `plot_altitude`, `plot_sky`, `plot_finder_image` (astroquery/SkyView), etc. matplotlib is imported lazily so the core package works without it.

## CI

- `.github/workflows/ci_tests.yml` calls the OpenAstronomy reusable tox workflow; each job is a tox env.
- `astropy/astroplan`'s branch protection on `main` requires CI checks **by job name** (`ci_tests / <name>`). Renaming or removing a job leaves PRs blocked until a repo admin updates the required checks, so call out any renames in the PR.
- Job names therefore describe roles, not versions ("oldest", "intermediate", "current" = `default_python`, and "latest" for devdeps, which moves to a new Python first). When bumping Python versions, change the tox envs and `default_python`, not the job names.

## Conventions

- Add an entry to `CHANGES.rst` under the unreleased version for user-facing changes, with the PR number in brackets (e.g. `[#605]`).
- Quantities with astropy units (`u.deg`, `u.minute`) and `astropy.time.Time` are used throughout; don't pass bare floats where the API expects Quantities.
