.. include:: references.txt

.. _astroplan:

**********************************
Observation Planning (`astroplan`)
**********************************

What is astroplan?
==================

**astroplan** is an open source Python package for observation planning for astronomers.

We're just getting started, this is not implemented yet, for now check out `Astropy`_,
`pyephem`_ or `skyfield`_.

The plan is to create a flexible toolbox for observation planning and scheduling both
for individual astronomers (professional or amateur) planning an observing night, as well
as observatories preparing long-term or nightly schedules, taking constraints (moon, slew
speeds, ...) into account.

We anticipate that astroplan will have the following required dependencies: 
* Numpy
* Astropy
* pytz

with potential optional dependencies including matplotlib and PyEphem or Skyfield.

A more detailed description of the scope and API proposal for ``astroplan`` will be posted in May/June 2015.

* Code, feature requests, bug reports, pull requests: https://github.com/astroplanners/astroplan
* Questions: http://groups.google.com/group/astropy
* Docs: https://astroplan.readthedocs.org/
* License: BSD-3
* We support Python 2.7 and 3.3+. (Python 2.6 and 3.2 or earlier are not supported.)

.. _astroplan_news:

News
====

* April 28, 2015: astroplan repo set up

.. _astroplan_docs:

General documentation
=====================

.. toctree::
  :maxdepth: 1

  getting_started/index

Reference/API
=============

.. automodapi:: astroplan
