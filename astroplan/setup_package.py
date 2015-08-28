# Licensed under a 3-clause BSD style license - see LICENSE.rst


def get_package_data():
    """Declare astroplan datafiles"""
    pdata = dict()
    pdata['astroplan.tests'] = ['coveragerc']
    pdata['astroplan.plots.tests'] = ['baseline_images/*.png']
    return pdata
