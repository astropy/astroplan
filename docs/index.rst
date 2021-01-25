.. include:: references.txt

.. _astroplan:

**********************************
Observation Planning (`astroplan`)
**********************************

What is astroplan?
==================

**astroplan** is an open source Python package to help astronomers plan
observations.

The goal of astroplan is to make a flexible toolbox for observation planning and
scheduling.  When complete, the goal is to be easy for Python beginners and new
observers to to pick up, but powerful enough for observatories preparing nightly
and long-term schedules.


Features:

* Calculate rise/set/meridian transit times, alt/az positions for targets at
  observatories anywhere on Earth
* Built-in plotting convenience functions for standard observation planning
  plots (airmass, parallactic angle, sky maps).
* Determining observability of sets of targets given an arbitrary set of
  constraints (i.e., altitude, airmass, moon separation/illumination, etc.).
* `Astropy <https://astropy.org>`__ powered!

Links
=====

* `Source code <https://github.com/astropy/astroplan>`_
* `Docs <https://astroplan.readthedocs.io/>`_
* `Issues <https://github.com/astropy/astroplan/issues>`_

License: BSD-3

.. _astroplan_docs:

General Documentation
=====================

.. toctree::
   :maxdepth: 2

   installation
   getting_started
   tutorials/index
   faq/index
   api
   changelog

.. _astroplan_authors:

Authors
=======

Maintainers
-----------
* `Jazmin Berlanga Medina, including contributions from Google Summer of Code 2015 <https://www.google-melange.com/gsoc/project/details/google/gsoc2015/jberlanga/5707702298738688>`_
* `Brett Morris, including contributions from Google Summer of Code 2015 <https://www.google-melange.com/gsoc/project/details/google/gsoc2015/bmmorris/5707702298738688>`_

Contributors
------------
* Christoph Deil
* Stephanie Douglas
* Eric Jeschke
* Adrian Price-Whelan
* Erik Tollerud
* Brigitta Sipocz
* Karl Vyhmeister
