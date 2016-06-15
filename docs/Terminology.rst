.. _Observation_Terminology:

***************
Terminology
***************

Scheduling
==================
Determined during the python in astronomy conference:

* observing block: the minimum schedulable block
* priority: number assigned to the observing block by the user (or TAC) that 
  defines its precedence within the set of blocks to be scheduled. Should probably
  also be on a 0->1 (0=no good, 1=good) scale, or rescaled within the scheduler
* constraint: sets limits and can yield a scores based on inputs. Scores can be 
  boolean or floats (0=no good, 1=good), with a flag in the Constraint call 
  that selects which is used
* score: the value returned by evaluating a constraint
* scorekeeper: assigns a cumulative score to the OB based on some function that 
  would be applied to all the individual scores (result should be in 0, 1)
* scheduler: the entity that consumes the results of the scorekeeper and the 
  observing blocks and produces a schedule
* Weights (“user defined”?): preferences for which constraints matter most 
  (i.e. I care more about getting dark time than getting a low airmass)
