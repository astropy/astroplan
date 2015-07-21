.. include:: references.txt

.. _astroplan:

**********************************
Observation Planning (`astroplan`)
**********************************

What is astroplan?
==================

**astroplan** is an open source Python package that helps astronomers plan
observations.

We're just getting started, so for now check out `Astropy`_, `PyEphem`_ or
`Skyfield`_.

**astroplan** is a flexible toolbox for observation planning and scheduling,
and is appropriate for astronomers at all skill levels.
It's also great for observatories preparing long-term or nightly schedules,
since **astroplan** can consider multiple constraints (i.e., airmass, moon,
slew speeds, etc.)

We anticipate that **astroplan** will have the following required dependencies:

* Python 2.7 or 3.3+ (Python 2.6 and 3.2 or earlier are not supported.)
* Numpy
* Astropy
* pytz

with potential optional dependencies including Matplotlib and PyEphem (or Skyfield).

**A first version is expected to roll out in late August 2015.**

Links
=====

* `Code, feature requests, bug reports, pull requests <https://github.com/astroplanners/astroplan>`_
* `Questions <http://groups.google.com/group/astropy>`_
* `Docs <https://astroplan.readthedocs.org/>`_

License: BSD-3

.. _astroplan_news:

News
====

* April 28, 2015: **astroplan** repo set up
* May/June 2015: API proposal posted on GitHub.

.. _astroplan_docs:

General documentation
=====================

.. toctree::
   :maxdepth: 1

   installation.rst
   getting_started.rst
   api.rst
