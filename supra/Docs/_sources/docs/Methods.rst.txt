.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

#######
Methods
#######

***********
Ray-Tracing
***********

******************************
Finding the Wave Release Point
******************************
To find the wave release point for a given station, many points are taken along the trajectory, and a ray-trace is calculated from every one of these points to the station. The initial direction of the acoustic path is then compared to the trajectory vector. Ideally, the initial direction of the path must be perpendicular to the trajectory vector in order for that ray to be the true acoustic path, but a tolerence is used in finding this.

.. note::
	If there are no winds or temperature gradients, the ray is a straight line from the trajectory to the station, and therefore is a lot easier to calculate where the ray is perpendicular to the trajectory vector. When winds are added, the acoustic path becomes curved, (and may not even arrive at the station), so the path is calculated, and if it hits the station (within a tolerence), the initial direction of this ray is used.

The angle tolerance may be adjusted to allow more ballistic arrivals at stations if none seem to be showing up. The minimum and maximum heights of wave release points may be adjusted to restrict to the known start and finish of the meteor.

.. warning::
	When the minimum height of the wave release points is set too low, the lower points along the meteor tend to release perpendicular to the trajectories, which may be unphysical. Tweaking these parameters per meteor is suggested.

************************
Slicing Atmospheric Data
************************

The path between the source and the detector are assumed to be linear for atmospheric purposes only. The linear path is divided up into N (default=10) points. Each point is moved to the closest 0.25x0.25 degree grid point. At each grid point, the Nth fraction of the atmospheric profile is taken from that point. The atmospheric profile fractions are then merged into the final atmospheric profile for the sound wave. The boundaries between atmospheric sections are interpolated in order to maintain a smooth atmospheric profile.

***********************
Finding the Supracenter
***********************

**********************
Finding the Trajectory
**********************

**********************
Perturbing the Weather
**********************

*****************
Finding the Yield
*****************



.. toctree::
   :maxdepth: 2

   





