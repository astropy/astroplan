"""
The scheduler is a very useful tool. It requires the input of observing blocks
and the time_range. Observing blocks consist of a target, a duration and constraints.
"global" constraints can also be defined and can significantly speed up the scheduler.
"""

# a few packages are needed to set up the proper objects
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astroplan
from astroplan.constraints import *
from astroplan import Observer, FixedTarget
from astroplan.scheduling import *

# create the time frame during which you observing
schedule_start = Time.now()+4*u.hour
schedule_end = schedule_start+24*u.hour

# define the observer
mro = Observer.at_site('mro')

# we also want to define how quickly the telescope can transition
# between different targets, filters, instruments etc.
# for using a telescope that slews at 1 degree/second,
# takes 2 minutes to change filters, and 10 to change instruments
slew_rate = 1*u.deg/u*second
transitioner = Transitioner(slew_rate, {filter: 2*u.minute, instrument: 10*u.minute})

# make the targets (this can also use FixedTarget(coords) if you make a
#                   `SkyCoord` object "coords")
targets = [FixedTarget.from_name('vega'),
           FixedTarget.from_name('Arcturus'),
           FixedTarget.from_name('Altair'),
           FixedTarget.from_name('Aldebaran'),
           FixedTarget.from_name('Sirius'),
           FixedTarget.from_name('Betelgeuse'),
           FixedTarget.from_name('Rigel'),
           FixedTarget.from_name('Castor')]
target_group = [FixedTarget.from_name('Alcor'),
                FixedTarget.from_name('M101')]

priorities = {target.name: i for i, target in enumerate(targets)}

# now `ObservingBlock`s and `BlockGroup`s need to be defined
visit_duration = 55*u.minute
blocks = []

# a `BlockGroup` contains one or more `ObservingBlocks` that are related

# create blocks for the targets  with the constraint that the airmass be less than 3
# for that target during its OB and a non-boolean score (since constraints default to a boolean)
# so that lower airmass gets a higher score. We also want the moon to be farther than 60 degrees
# away, but lower airmass is better than the moon being farther, so we give them weights.
for target in targets:
    block = BlockGroup(target, visit_duration, priority=priorities[target.name],
                       constraints=[AirmassConstraint(max=3, boolean_score=False),
                                    MoonSeparationConstraint(60*u.deg, boolean_score=False)],
                       constraint_weights=[1, .2],
                       instrument_configurations=[{filter: 'v'}]
                       )
    blocks.append(block)

# create a `BlockGroup` to make sure Alcor and M101 are observed consecutively
# we don't care how good the airmass is for them, just that it is below 3 and
# that the observations are done during grey time.
# we also want the group prioritized after the single targets and
# to observe M101 immediately before Alcor (i.e. order = [1, 0]
blocks.append(BlockGroup(target_group, [visit_duration for target in target_group],
                         priority=max(priorities)+1, order=[1, 0],
                         constraints=[AirmassConstraint(3), MoonIlluminationConstraint.grey()],
                         instrument_configurations=[{filter: 'b'}, {filter: 'v'}],
                         transitioner = transitioner
                         ))

# The observatory is optical so all ObservationBlocks need to be scheduled at night
# it also can only point between 10 and 80 degrees altitude.
observer_constraints = [AtNightConstraint(), AltitudeConstraint(10*u.deg, 80*u.deg)]

# now create a scheduler object, in this case we want the scheduler which uses
# the Greedy method for scheduling
scheduler = GreedyScheduler(schedule_start, schedule_end, constraints=observer_constraints,
                            observer=mro, transitioner=transitioner)

# run the scheduler on the defined blocks
schedule = scheduler(blocks)

# to see the chronological list of the OBs and when they are scheduled for
print(schedule.scheduled)

# to see the above, but in table form instead (with columns of target, start, end, target_ra, target_dec, score, etc.)
print(schedule.to_table)

# plotting when the targets were scheduled and the airmass of the targets.
plot_schedule_airmass(schedule)

# to see how much of the time is spent observing
print(schedule.up_time, schedule.down_time)
