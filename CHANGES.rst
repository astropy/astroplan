0.11 (unreleased)
-----------------

- Fix ``LocalTimeConstraint`` comparing its limits against UTC instead of the
  observer's local time (or the timezone attached to ``min``/``max``). Results
  now change for any observer not in UTC. The constraint also now works on
  multidimensional time grids and no longer overwrites ``min`` when only
  ``max`` is given. [#466]

- Replace ``pytz`` with the standard library's ``zoneinfo``. ``pytz`` is no
  longer a dependency; ``tzdata`` is required instead so timezone data is
  available on all platforms. ``Observer`` now creates ``zoneinfo.ZoneInfo``
  objects when ``timezone`` is given as a string, so ``Observer.timezone`` no
  longer has pytz-specific methods such as ``localize``. Timezones created with
  ``pytz`` are still accepted. [#601]

- Minimum Python version is now 3.11. Also bumped minversion
  of dependencies. [#605]

- Minimum Python version is now 3.12, and Python 3.14 is now tested. Minimum
  versions of numpy (1.26.0) and matplotlib (3.7.3) were bumped to the oldest
  releases available for Python 3.12.

- The ``oldestdeps`` test environment now uses ``tox-uv`` with uv's ``lowest``
  resolution strategy, so it tests the lower bounds declared in
  ``pyproject.toml`` (plus lower bounds for transitive dependencies in
  ``oldestdeps-constraints.txt``) instead of hand-maintained pins. This
  raised the minimum version of pytest (7.4.0), and
  added minimum versions for astroquery (0.4.7), pytest-astropy (0.11.0) and
  pytest-mpl (0.16.1). [#631]

- Add a ``dev`` extra that installs the ``all`` and ``test`` extras plus
  ``tox`` and ``tox-uv``, for running the test environments locally. [#641]

- Fix CI: bump the pinned ``OpenAstronomy/github-actions-workflows`` reusable
  workflow to v3.0.3 (fixes a Python-3.12-patch-version matching bug that was
  intermittently breaking the "Load tox environments" and wheel-building
  steps), ensure the coverage data file lands where the new Codecov upload
  step expects it, fix stale ``matplotlib.pyplot.plot_date`` docstring
  references left over from the switch to ``matplotlib.pyplot.plot``, fix the
  now fully-removed ``Axes.plot_date`` call in the IERS FAQ's example plot,
  fix a dead documentation link, and ignore ADS's bot-blocking 405 responses
  in the linkcheck.

0.10.1 (2024-08-13)
-------------------

- Compatibility with numpy 2.0. [#587]

0.10 (2024-04-04)
-----------------

- Fix compatibility with astropy v6.0.

0.9.1 (2023-09-20)
------------------

- Fix bug when ``FixedTarget`` objects are passed to methods that calculate
  lunar coordinates. [#568]

0.9 (2023-07-27)
----------------

- Fix time range in ``months_observable`` to not be only in 2014. Function now
  accepts argument ``time_range`` and defaults to the current year. [#458]

- Fix ``Observer`` not having longtitude, latitude, and elevation parameters
  as class attributes. They are now properties calculated from the ``location``.

- Documentation revisions and theme update [#563]

0.8 (2021-01-26)
----------------

- Fix Read The Docs compatibility [#497]

- Move to APE 17 infrastructure, change to github actions [#493]

- Update conda channel in favor of conda-forge [#491]

- Fix for astropy cache compatibility [#481]

0.7 (2020-10-27)
----------------

- Fix compatibility with Astropy 4.X


0.6 (2019-10-08)
----------------

- Added documentation for reproducing MMTO sun rise/set times [#434]

- Deprecation of ``MAGIC_TIME`` variable, which used to be returned for targets
  that don't rise or set [#435]

- Replace deprecated astroquery service [#431]

- Fix for the broken IERS patch [#418, #425]

- Add ``GalacticLatitudeConstraint`` to constrain the galactic latitudes of
  targets. This can be useful for planning surveys for which crowding due to
  Galactic point sources is an issue. [#413]


- Add ``n_grid_points`` keyword argument to rise/set/transit functions which
  allows usersto trade off precision for speed. [#424]

0.5 (2019-07-08)
----------------

- ``observability_table`` now accepts scalars as ``time_range`` arguments, and
  gives ``'time observable'`` in this case in the resulting table. [#350]

- Bug fixes [#414, #412, #407, #401]

0.4 (2017-10-23)
----------------

- Added new ``eclipsing`` module for eclipsing binaries and transiting
  exoplanets [#315]

- Fixes for compatibility with astropy Quantity object updates [#336]

- Better PEP8 compatibility [#335]

- Using travis build stages [#330]

0.3 (2017-09-02)
----------------

- ``Observer.altaz`` and ``Constraint.__call__`` no longer returns an (MxN) grid
  of results when called with M ``target``s and N ``times``. Instead, we attempt
  to broadcast the time and target shapes, and an error is raised if this is not
  possible. This change breaks backwards compatibility but an optional argument
  ``grid_times_targets`` has been added to these methods. If set to True,
  the old behaviour is recovered. All ``Observer`` methods for which it is
  relevant have this optional argument.

- Updates for compatibility with astropy v2.0 coordinates implementation
  [#311], updates to astropy-helpers [#309], fix pytest version [#312]

0.2.1 (2016-04-27)
------------------

- Internal changes to the way calculations are done means that astropy>=1.3 is required [#285]

- Fixed bug when scheduling block list is empty [#298]

- Fixed bug in Transitioner object when no transition needed [#295]

- Update to astropy-helpers 1.3.1 [#294] and compatibility fixes for astropy 1.3 [#283]


0.2 (2016-09-20)
----------------

- Fixed bug arising from changes to distutils.ConfigParser [#177, #187, #191]

- Removed the sites module from astroplan, since it was ported to astropy [#168]

- Removed dependence on PyEphem, now using jplephem for the solar system
  ephemeris [#167]

- New API for scheduling observations (still in development)

- New ``plot_finder_image`` function makes quick finder charts [#115]

- Updates to astropy helpers and the package template [#177, #180]
