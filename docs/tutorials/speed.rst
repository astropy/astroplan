.. _speed:

*********************
Speeding up astroplan
*********************

Some users of astroplan may find it useful to trade-off a bit of precision
in the rise/set/transit times of targets in exchange for computational
efficiency. In this short tutorial, we show you how to speed up astroplan
in exchange for a bit of time precision, which is especially useful when
planning many observations over a long period of time.

Rise/set/transit times
======================

The rise, set, and transit time methods on the `~astroplan.Observer` object
take an optional keyword argument called ``N`` as of astroplan version 0.6
(in earlier versions of astroplan, ``N`` is fixed to 150). To understand ``N``
you first need to know how target rise/set/transit times
are computed in astroplan.

Astroplan computes rise/set times relative to a given reference time by
computing the altitude of the target on a grid which spans a period of 24 hours
before/after the reference time. The grid is then searched for
horizon-crossings, and astroplan interpolates between the two nearest-to-zero
altitudes to approximate the target rise/set times.

The ``N`` keyword argument dictates the number of grid points on which to
compute the target altitude. The larger the ``N``, the more precise the
rise/set/transit time will be, but the operation also becomes more
computationally expensive. As a general rule of thumb, if you choose ``N=150``
your rise/set time precisions will be precise to better than one minute; this
is the default if you don't specify ``N``. If you choose ``N=10`` you'll get
significantly faster rise/set time computations, but your precision degrades to
better than five minutes.

Examples
========

Let's see some simple examples. We can compute a very accurate rise time for
Sirius over Apache Point Observatory, by specifying ``N=1000``:

.. code-block:: python

    >>> from astroplan import Observer, FixedTarget
    >>> from astropy.time import Time

    >>> time = Time('2019-01-01 00:00')
    >>> sirius = FixedTarget.from_name('Sirius')
    >>> apo = Observer.at_site('APO')

    >>> rise_time_accurate = apo.target_rise_time(time, sirius, N=1000)
    >>> rise_time_accurate.iso  # doctest: +SKIP
    '2019-01-01 01:52:13.393'

That's the rise time computed on a grid of 1000 altitudes in a 24 hour period,
so it should be very accurate. Now let's compute a lower precision, but much
faster rise time, using ``N=10`` this time:

.. code-block::

    >>> rise_time_fast = apo.target_rise_time(time, sirius, N=10)
    >>> rise_time_fast.iso  # doctest: +SKIP
    '2019-01-01 01:54:09.946'

You can see that the rise time returned by
`~astroplan.Observer.target_rise_time` with ``N=10`` is only two minutes
different from the prediction with ``N=1000``, so it looks like we haven't lost
much precision despite the drastically different number of grid points.