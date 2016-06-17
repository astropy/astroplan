.. _Observation_Terminology:

***********
Terminology
***********

Scheduling
==========

Scheduling observations is a common task in astronomy, but there is a lack of
formal terminology. Here we define the terminology as it will be used in astroplan
for scheduling tasks. This is for consistency within the code and the documentation.

The terms used are as follows:

* **visit**: a period of time spent observing a single target
* **observing block** (OB): the minimum schedulable block (tentative: a group of visits)
* **priority**: number assigned to the **observing block** by the user that 
  defines its precedence within the set of blocks to be scheduled. Should probably
  also be on a 0->1 (0=no good, 1=good) scale, or rescaled within the scheduler
* **rank**: a **priority** defined by a TAC rather than a user
* **constraint**: sets limits and can yield **scores** based on inputs. Currently
  only boolean scores are implemented in the code.
* **score**: the value returned by evaluating a **constraint**. can be 
  boolean or floats between 0 and 1 inclusive (0=no good, 1=good), with a flag in the 
  ``Constraint`` call that selects which is used
* **scorekeeper**: assigns a cumulative score to the OB based on some function that 
  would be applied to all the individual **scores** (result should be in 0, 1)
* **scheduler**: the entity that consumes the results of the **scorekeeper** and the 
  **observing blocks** and produces a schedule
* **Weights** (“user defined”?): preferences for which constraint's **scores** matter most 
  (i.e. I care more about getting dark time than getting a low airmass)


This is a list of possible schedulers, none are implemented yet. Once implemented they
will have different methods for creating their schedule. Below is a list of ideas for
the scheduler descriptors and the related scheduler would work. Each scheduler will be
able to be called with ``DescriptorScheduler(args)`` using the descriptor defined below.

* **Sequential**: starts from the beginning of the time range and schedules the OB
  with the best **score**, that hasn't already been scheduled, for that time.
* **Priority**: starts from the highest **priority** OB and schedules it at the time
  where it has its highest **score**. Then schedules the next-highest priority without
  overlapping