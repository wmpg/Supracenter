.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#####
Tools
#####

Trajectory Interpolator
=======================

The Trajectory Interpolator will produce a .csv of 4-Dimensional points of the trajectory defined in the current event .bam file. Note that currently the points produced are only an approximation. The .csv produced by this tool may be directly imported into the event .bam file to quickly plot a bunch of test fragmentations along the trajectory. It is intended to use this along with the Fragmentation Staff to identify the heights of possible fragmentation arrivals for a given station.

Stratospheric Circulation Index
===============================

Used to estimate wind effects over long ranges. To use, go to the "Fetch Atmosphere" tab, and select "SCI Wind Index". The terminal will show the average U- and V- component of the wind values between 40 and 60 km of the currently loaded wind data (if winds are off, then HWM14, if winds are on, then ECMWF). The SCI is calculated as follows: First the wind data is interpolated from the beginning to end points as a straight line, taking data from the nearest 0.5 degree grid point. This straight line is merged into a single sounding profile with height, which roughly matches the path selected. The data is then cubic splined to create a smooth profile. The mean of the values between 40 and 60 km are taken to get the U- and V- component of the SCI. The angle from the beginning to end point (taking into account curvature of the Earth), is taken, and used to calculate the SCI component in the direction of travel.

.. toctree::
   :maxdepth: 2

   





