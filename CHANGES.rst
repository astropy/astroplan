0.3 (unreleased)
----------------

- Internal changes to the way calculations are done means that astropy>=1.3 is required [#285]

0.2.1 (2016-04-27)
------------------

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
