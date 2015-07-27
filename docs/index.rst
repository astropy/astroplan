.. include:: references.txt

.. _astroplan:

**********************************
Observation Planning (`astroplan`)
**********************************

What is astroplan?
==================

**astroplan** is an open source Python package that helps astronomers plan
observations.

We're still in development, so for now check out `Astropy`_, `PyEphem`_ or
`Skyfield`_.

Our aim is to make astroplan a flexible toolbox for observation planning and
scheduling.  We want astroplan to be easy for Python and observing beginners
to pick up, but powerful enough for observatories preparing nightly and
long-term schedules.

Features:

* `Astropy`_ powered!
* Ability to take multiple constraints into account (i.e., airmass, moon, slew speeds, etc.) to determine visibility of objects.
* Built-in plotting functions (airmass, parallactic angle, sky maps)

We anticipate that **astroplan** will have the following required dependencies:

* Python 2.7 or 3.3+ (Python 2.6 and 3.2 or earlier are not supported.)
* Numpy
* Astropy
* pytz

with potential optional dependencies including Matplotlib and PyEphem
(or Skyfield).

**A first version is expected to roll out in late August 2015.**

Links
=====

* `Code, feature requests, bug reports, pull requests <https://github.com/astropy/astroplan>`_
* `Questions <http://groups.google.com/group/astropy>`_
* `Docs <https://astroplan.readthedocs.org/>`_

License: BSD-3

.. _astroplan_news:

News
====

* ???: Astroplan **0.1** release
* May - August 2015: Initial package implemented by Jazmin Berlanga and Brett Morris
  in `GSoC 2015 <https://www.google-melange.com/gsoc/homepage/google/gsoc2015>`.

.. _astroplan_docs:

General Documentation
=====================

.. toctree::
   :maxdepth: 1

   installation.rst
   getting_started.rst
   tutorials.rst
   howto.rst
   api.rst
