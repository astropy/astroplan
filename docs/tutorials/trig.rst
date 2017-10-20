.. _trig_approx_tutorial:

*****************************************
Fast trigonometric altitude approximation
*****************************************

Contents
========

* :ref:`trig-background`
* :ref:`trig-trig_approx`

.. _trig-background:

Background
==========

astroplan ordinarily computes the rise and set times of an object by
transforming the sky coordinate of the object (e.g.~ICRS, galactic, etc.) into a
grid of altitude-azimuth coordinates for that target as seen by an observer at a
specific location on the Earth, at 10 minute intervals over a 24 hour period.
The rise or set time is then computed by linear interpolation between the two
coordinates nearest to zero. The meridian/anti-meridian transit time is computed
similarly; it takes a numerical derivative of the altitudes before searching for
the appropriate zero crossing.

We chose to compute rise and set times with a grid-search to maximize accuracy,
rather than speed. In particular, we sought to preserve the astropy
altitude-azimuth coordinate transformation which accounts for atmospheric
refraction. If you have a need for speed, however, there's an alternative...

.. _trig-trig_approx:

The ``trig_approx`` option
==========================

A simple trigonometric algorithm is included within the `~astroplan.Observer`
object for computing rise and set times without accounting for refraction, which
yields a 3x speed gain in exchange for diminished accuracy. The fast,
trigonometric altitude approximation yields rise/set time precisions <10 min
when compared against equivalent computations in PyEphem.

One contribution to the difference in rise/set times between the packages is the
different interpretations of time systems in astroplan and pyephem. astroplan
will use the fast trigonometric altitude approximation when an
`~astroplan.Observer` object is initialized with the ``trig_approx`` option:


.. code-block:: python

    >>> from astroplan import Observer
    >>> keck_fast = Observer.at_site("Keck", trig_approx=True)

When the ``trig_approx == True``, astroplan will use the fast trigonometric
approximation when the following methods are called:

- `~astroplan.Observer.target_rise_time`
- `~astroplan.Observer.target_set_time`
- `~astroplan.Observer.target_meridian_transit_time`
- `~astroplan.Observer.target_meridian_antitransit_time`
- `~astroplan.Observer.target_is_up`
- `~astroplan.Observer.is_night`


.. code-block:: python

    >>> from astropy.time import Time
    >>> reference_time = Time('2015-06-07 00:00')

    >>> from astroplan import FixedTarget
    >>> target = FixedTarget.from_name('Sirius')

    >>> approx_rise = keck_fast.target_rise_time(reference_time, target)
    >>> print(approx_rise.iso) # doctest: +FLOAT_CMP
    2015-06-06 18:32:31.960

This rise time result will be different from the full-precision calculation by
a small amount:

    >>> keck_precise = Observer.at_site('Keck')  # by default, trig_approx=False
    >>> precise_rise = keck_precise.target_rise_time(reference_time, target)
    >>> print(precise_rise.iso) # doctest: +SKIP
    2015-06-06 18:33:14.136

We can check here that the difference in times is small:

    >>> import astropy.units as u
    >>> abs(precise_rise - approx_rise) < 5*u.min
    True
