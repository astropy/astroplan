0.11 (unreleased)
-----------------

- Fix ``Observer.moon_phase`` failing when called without a time. It now
  uses the current time, like ``Observer.tonight``. ``Observer.moon_illumination``
  now also defaults to the current time, and ``Observer.tonight`` accepts
  anything ``Time`` accepts (e.g. a string), like the other ``Observer``
  methods. [#666]

- Add ``HourAngleConstraint`` to constrain a target's hour angle, e.g. to stay
  within the hour angle limits of an equatorial mount. ``min`` and ``max`` are
  angles; if ``min`` is greater than ``max``, the window wraps through 12h. [#560]

- Fix rise, set, transit, twilight and midnight times coming back masked when
  the next or previous event is slightly more than a day away (e.g. the Sun
  when day length is changing, or the Moon's transit) or falls exactly at the
  reference time. The search now covers one day plus a margin. Also,
  ``which='nearest'`` now returns the event on the other side of the reference
  time instead of a masked time when only one side has an event, and only
  warns when neither side does. [#656]

- Fix rise, set and meridian transit times with ``which='next'`` (or
  ``'nearest'``) skipping an event shortly after the reference time and
  returning the same event one sidereal day later instead. [#651]

- Add ``Observer.get_body``, a cached wrapper around
  ``astropy.coordinates.get_body`` for the observer's location, and use it for
  the Sun and Moon positions in ``Observer`` and in the Sun/Moon separation and
  ``AtNightConstraint`` constraints. ``moon_illumination`` and
  ``moon_phase_angle`` accept an ``observer`` argument to reuse its cached
  positions (shifted back to the geocentre, so results are unchanged);
  ``Observer.moon_illumination``, ``Observer.moon_phase`` and
  ``MoonIlluminationConstraint`` use it [#662]. The Sun
  helpers (``sun_rise_time``, ``sun_set_time``, ``noon``, ``midnight``,
  ``sun_altaz``, ``is_night`` and ``altaz`` with the Sun) now start from
  ``get_body('sun', ...)`` instead of ``get_sun``; because the AltAz transform
  already applied the topocentric correction, computed altitudes change by
  less than 0.03 arcsec and rise/set times by about a millisecond. [#578]

- Add ``MeridianSeparationConstraint`` to constrain a target's angular
  separation from the meridian, e.g. to avoid pier flips on German
  equatorial mounts. [#589]

- ``SunSeparationConstraint`` and ``MoonSeparationConstraint`` no longer emit
  ``NonRotationTransformationWarning``. The minimum required astropy version
  is now 6.1. [#649]

- Add ``hours_value`` option to ``plot_sky`` to annotate each plotted
  position with a label (e.g. the hour). A per-point ``c`` array in
  ``style_kwargs`` is now masked along with below-horizon positions so
  colormaps stay matched to the plotted points. [#617]

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

- Fix ``PriorityScheduler`` raising ``block already occupied`` for observing
  blocks with durations at or below ``time_resolution``. Block durations are
  now rounded up to the scheduler's time grid; for example, a 55-minute block
  with 2-minute resolution now lasts 56 minutes instead of 54. [#550]

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

- ``Observer.parallactic_angle`` and ``Observer.target_hour_angle`` now
  transform the target to the apparent topocentric ``HADec`` frame before
  computing the hour angle, instead of subtracting the target's input-frame
  (e.g. ICRS) right ascension from the local sidereal time. This accounts for
  precession, nutation, aberration and topocentric parallax, so results
  change by ~0.1-0.4 deg for typical targets at current epochs, and by more
  close to the pole or zenith. ``parallactic_angle`` now also accepts
  ``MoonFlag`` and ``SunFlag`` targets. Its ``kind`` and ``model`` arguments
  are deprecated and ignored. [#644]

- Fix rise/set times always being returned as masked ``Time`` objects, even
  when no times are masked. Newer astropy keeps such ``Time`` objects masked,
  which broke coordinate transforms that use them (e.g. ``Observer.tonight``).
  [#643]

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

- Fix ``PhaseConstraint`` raising a ``TypeError`` when ``min`` or ``max`` is
  not given, instead of defaulting to the full phase range. [#639]

- Add ``AltAzTarget`` for targets defined by a fixed altitude and azimuth
  at a given observing location, with support for time-dependent targets
  in calculations, constraints and scheduling. ``Schedule.to_table()``
  now also includes target type and target information. [#634]

- With ``grid_times_targets=True``, a single time and multiple targets
  now return an (M, 1) array instead of (M,), including when passing a
  multi-coordinate ``SkyCoord``. [#634]

- ``PriorityScheduler`` now rounds constraint scores to six decimal places
  when choosing observation times, making scheduling more consistent
  and preferring earlier slots when scores are nearly equal. [#634]

- ``Scorer.targets`` now stores the original targets as a list instead
  of combining them into a ``SkyCoord``. [#634]

- Add ``SGP4SatelliteTarget`` for satellites propagated from TLE or OMM
  elements using the optional ``sgp4`` dependency. Satellite targets can be
  used with observer calculations, constraints and scheduling, and
  ``Schedule.to_table()`` displays their catalog number and element epoch.
  [#566]

- Improve performance of night-time shading in ``plot_airmass``,
  ``plot_altitude``, and ``plot_schedule_airmass`` by sampling the Sun's
  altitude and interpolating twilight boundaries. Handle missing
  sunrise and twilight events at high latitudes, fix shading across
  multiple nights, and shade empty schedules. The ``show_night`` option
  now uses the darker, graded twilight shading of ``plot_airmass``.
  [#559]

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
