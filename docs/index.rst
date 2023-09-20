.. include:: references.txt

.. _astroplan:

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
* `astropy <https://astropy.org>`__ powered!

Links
=====

* `Source code <https://github.com/astropy/astroplan>`_
* `Docs <https://astroplan.readthedocs.io/>`_
* `Issues <https://github.com/astropy/astroplan/issues>`_

License: BSD-3

.. _astroplan_docs:

Documentation
=============

.. toctree::
   :maxdepth: 2

   installation
   getting_started
   tutorials/index
   faq/index
   api
   changelog

Maintainers
+++++++++++
* `Brett Morris, including contributions from Google Summer of Code 2015 <https://www.google-melange.com/gsoc/project/details/google/gsoc2015/bmmorris/5707702298738688>`_

Attribution
+++++++++++

If you use astroplan in your work, please cite `Morris et al. 2018 <https://ui.adsabs.harvard.edu/abs/2018AJ....155..128M/abstract>`_:

.. code :: bibtex

  @ARTICLE{2018AJ....155..128M,
         author = {{Morris}, Brett M. and {Tollerud}, Erik and {Sip{\H{o}}cz}, Brigitta and {Deil}, Christoph and {Douglas}, Stephanie T. and {Berlanga Medina}, Jazmin and {Vyhmeister}, Karl and {Smith}, Toby R. and {Littlefair}, Stuart and {Price-Whelan}, Adrian M. and {Gee}, Wilfred T. and {Jeschke}, Eric},
          title = "{astroplan: An Open Source Observation Planning Package in Python}",
        journal = {\aj},
       keywords = {methods: numerical, methods: observational, Astrophysics - Instrumentation and Methods for Astrophysics},
           year = 2018,
          month = mar,
         volume = {155},
         number = {3},
            eid = {128},
          pages = {128},
            doi = {10.3847/1538-3881/aaa47e},
  archivePrefix = {arXiv},
         eprint = {1712.09631},
   primaryClass = {astro-ph.IM},
         adsurl = {https://ui.adsabs.harvard.edu/abs/2018AJ....155..128M},
        adsnote = {Provided by the SAO/NASA Astrophysics Data System}
  }