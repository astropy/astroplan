"""
The scheduler is a very useful tool. It requires the input of observing blocks
and the time_range. Observing blocks consist of a target, a duration and constraints.
"global" constraints can also be defined and can significantly speed up the scheduler.
"""

#a few packages are needed to set up the proper objects
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astroplan
from astroplan.constraints import *
from astroplan import Observer, FixedTarget
from astroplan.scheduling import *

#create the time frame during which you observing
start_time = Time.now()+4*u.hour
end_time = start_time+24*u.hour

#define the observer
mro = Observer.at_site('mro')

#make the targets (this can also use FixedTarget(coords) if you make a
#                   SkyCoord object coords)
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
#now ObservingBlocks and ScheduleBlocks need to be defined
#priorities = {target.name : i for i, target in enumerate(targets)}#optional code
visit_duration = 55*u.minute
blocks = []
#create a block for each target with the constraint that the airmass be less than 3
#for that target during its OB and a non-boolean score so that lower airmass gets a
#higher score.
for target in targets:
    blocks.append(ObservingBlock(target,visit_duration
                                 #,priority = priorities[target.name],
                                 constraints = [AirmassConstrant(3,bool_score = False)]
                                 ))
#create a ScheduleBlock to make sure Alcor and M101 are observed consecutively
#we don't care how good the airmass is for them, just that it is below 3 and 
#that the observations are done during grey time.
blocks.append(ScheduleBlock(target_group,
                            durations = [54*u.minute for target in target_group],
                            #priority = None,
                            #order = [1,0]
                            constraints = [AirmassConstraint(3),
                                           MoonIlluminationConstraint.grey()]
                            ))

#The observatory is visual so all ObservationBlocks need to end up at night
#it also can only point between 10 and 80 degrees altitude. (some of this could
#end up as parts of the Observer object)
constraints = [AtNightconstraint(),
               AltitudeConstraint(10*u.deg,80*u.deg)]

#now create a scheduler object, in this case we want the scheduler which uses
#the Greedy method for scheduling
scheduler = GreedyScheduler(start_time, end_time, constraints=constraints, observer = mro)

#run the scheduler on the defined blocks
schedule  = scheduler(blocks)

#to see the chronological list of the OBs and when they are scheduled for
print(schedule.scheduled)
#to see the entire schedule including transitions and open slots
print(schedule)
#to see how much of the time is spent observing
print(schedule.up_time, schedule.down_time)
