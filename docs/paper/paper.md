title: 'astroplan: an astropy affiliated package for observation planning'
tags:
  - astronomy
  - Python
  - astropy-affiliated
authors:
 - name: Brett M. Morris
   orcid: 0000-0003-2528-3409
   affiliation: 1
 - name: Karl Vyhmeister
   orcid: 
   affiliation: 2
 - name: Stuart Littlefair
   orcid: 0000-0001-7221-855X
   affiliation: 3
 - name: Brigitta Sipocz
   orcid: 
   affiliation: 4
 - name: Erik Tollerud
   orcid: 0000-0002-9599-310X
   affiliation: 5
 - name: Jazmin Berlanga Medina
   orcid: 
   affiliation: 6
 - name: Adrian Price-Whelan
   orcid: 0000-0003-0872-7098
   affiliation: 7
 - name: Wilfred Tyler Gee
   orcid: 0000-0002-2931-7605
   affiliation: 8
 - name: Stephanie Douglas
   orcid: 0000-0001-7371-2832
   affiliation: 9
 - name: Aniket Pandey
   orcid: 
   affiliation: 10
 - name: Kelle Cruz
   orcid: 0000-0002-1821-0650
   affiliation: 11
 - name: Larry Bradley
   orcid: 0000-0002-7908-9284
   affiliation: 5
 - name: Thomas Robitaille
   orcid: 0000-0002-8642-1329
   affiliation: 12
 - name: Eric Jeschke 
   orcid: 
   affiliation: 13
affiliations:
 - name: University of Washington
   index: 1
 - name: California Institute of Technology
   index: 2
 - name: University of Sheffield
   index: 3
 - name: University of Hertfordshire
   index: 4
 - name: Space Telescope Science Institute
   index: 5
 - name: Purdue University
   index: 6
 - name: Princeton University
   index: 7
 - name: University of Hawaii at Hilo
   index: 8
 - name: Harvard-Smithsonian Center for Astrophysics
   index: 9
 - name: Indian Institute of Technology Kanpur
   index: 10
 - name: Hunter College
   index: 11
 - name: (freelance)
   index: 12
 - name: National Astronomical Observatory of Japan 
   index: 13

date: 2 September 2017
bibliography: paper.bib

---

# Summary

Astropy is the core Python package for fundamental astronomical computations. An 
ecosystem of astropy affiliated packages provides efficient access to more 
specific astronomical tasks. We present astroplan – an astropy affiliated 
Python package for observation planning and scheduling. 

The astroplan API is designed to provide efficient access to common 
observational quantities without requiring a detailed understanding of 
astropy's implementation of coordinate systems. Common astroplan calculations 
include celestial rise, set, and meridian transit times, and altitude-azimuth 
coordinates for targets at observatories anywhere on Earth. We also provide 
convenience functions for common observational plots, such as airmass and 
parallactic angle plots, and basic sky maps.

Astroplan also defines a generic framework for observation scheduling. Users 
can determine whether or not a target is observable given a variety of 
observing "constraints", such as airmass limits, time ranges, moon 
illumination/separation ranges, and more. A few schedulers are included which
divide observing time among a list of targets, given observing constraints on 
those targets.

# References