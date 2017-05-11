0.3 (unreleased)
----------------

- No changes yet

- ``Observer.altaz`` and ``Constraint.__call__`` no longer returns an (MxN) grid of results when
  called with M ``target``s and N ``times``. Instead, we attempt to broadcast the time and target
  shapes, and an error is raised if this is not possible. This change breaks backwards compatibility
  but an optional argument ``grid_times_targets`` has been added to these methods. If set to true,
  the old behaviour is recovered.

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
