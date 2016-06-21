from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.utils import minversion

from ..observer import Observer
from ..target import FixedTarget
from ..constraints import (AirmassConstraint)
from ..scheduling import (ObservingBlock,  PriorityScheduler)

APY_LT104 = not minversion('astropy', '1.0.4')

vega = FixedTarget(coord=SkyCoord(ra=279.23473479*u.deg, dec=38.78368896*u.deg),
                   name="Vega")
rigel = FixedTarget(coord=SkyCoord(ra=78.63446707*u.deg, dec=8.20163837*u.deg),
                    name="Rigel")
polaris = FixedTarget(coord=SkyCoord(ra=37.95456067*u.deg,
                                     dec=89.26410897*u.deg), name="Polaris")


def test_scheduler_hashable():
    # if it runs into an unhashable, this should error
    start_time = Time.now()
    end_time = start_time+48*u.hour
    mdm = Observer.at_site('mdm')
    targets = [vega, rigel, polaris]
    blocks = []
    for i, t in enumerate(targets):
        # i is used as the priority
        blocks.append(ObservingBlock(t, 55*u.min, i))
    constraints = [AirmassConstraint(3, boolean_constraint=False)]
    scheduler = PriorityScheduler(start_time, end_time,
                                  constraints=constraints, observer=mdm)
    schedule = scheduler(blocks)
    pass
