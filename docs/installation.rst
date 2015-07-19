.. include:: references.txt

.. _installation:

************
Installation
************

Requirements
============

**astroplan** requires Python 2.7 or 3.3+ (2.6 and 3.2 or earlier are not
supported) as well as the following packages:

* `Numpy`_
* `Matplotlib`_
* `Astropy`_
* `pytz`_

First-time Python users may want to consider an all-in-one Python installation
package, such as the `Anaconda Python Distribution
<http://continuum.io/downloads>`_ which provides all of the above dependencies.

Optional packages:

* `PyEphem`_
* `Skyfield`_

Installation
============

You can install the stable version of astroplan with::

    pip install astroplan

Alternatively, you can install the latest developer version of astroplan by
cloning the git repository::

    git clone https://github.com/astroplanners/astroplan

...then installing the package with::

    cd astroplan
    python setup.py install

Testing
=======

If you want to check that all the tests are running correctly with your Python
configuration, start up python, and type::

    import astroplan
    astroplan.test()

If there are no errors, you are good to go!

More
====

astroplan follows `Astropy`_'s guidelines for affiliated packages--installation
and testing for the two are quite similar! Please see Astropy's
`installation page <http://astropy.readthedocs.org/en/latest/install.html>`_
for more information.
