0.8 (unreleased)
----------------

- No changes yet.

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
