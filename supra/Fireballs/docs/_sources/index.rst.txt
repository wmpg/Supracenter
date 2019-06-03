.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.



Western Meteor Physics Group Seismic and Infrasonic Meteor Program User Manual
******************************************************************************

.. figure::  wmpl.png
   :align:   center




---------------------
Supracenter Manual Search
---------------------
Chose the 4-D fragmentation point, and input the station data and picks, created by the Make Picks window, and hit search. The program will return the error in that point, as well as the residuals to each station, based off of the difference in arrival times from the manual Supracenter to the expected pick arrival time.

---------------------
Supracenter PSO Search
---------------------
Chose a 4-D bounding box to search in, giving the 3-D spatial barrier and a 1-D temporal barrier. Give the station picks file, created by the Make Picks window. The program will search for the optimal Supracenter within the 3-D bounding box, and try to optimize the for the best solution with the lowest residuals. The residuals are calculated by the difference between the expected arrival time of a trial point, and the arrival time given in the picks file.

---------------------
Seismic Trajectory
---------------------
Functional, but not thouroughly tested, use at your own risk!

====================
Files
====================

--------------------
INI Files
--------------------

ini files contain the information for a given fireball event. These are loaded into each part of the program with the parameters of the fireball. For most parts of the program, not all of the variables in the ini files need to be set. For example, in the Supracenter solver, the ballistic search parameters do not need to be filled out.

If the terminal returns an "INI ERROR" while running the code, try setting [General] debug = True, which will print out what variables are being detected by the program.

--------------------
Station Picks
--------------------


--------------------
Atmosphere Files
--------------------


--------------------
All Times Files
--------------------

These files contain the estimated arrival times for each perturbation at each station for all fragmentation and ballistic arrivals. This .npy file is used to quickly save and load data and prevent long loading times. The data can be called from a numpy array with the data loaded in as follows:

[perturbation #, station #, 0 - ballistic/ 1 - fragmentation, frag number (0 if using ballistic)]

The intended use of this is to save estimates for arrival times that take a long time to load, such as perturbation calculations. It is recommended that each .npy file which contains perturbations is saved in a proper place so that it may be reused.



====================
Errors
====================

There definitely are a bunch, but this section is a WIP

====================
Extra Stations
====================

To add a station manually, the station name, location and mseed file must be input into the .ini file. Most seismic waveform databases online will allow for an mseed file to be downloaded. The file must be put in with the rest of the mseed files, which is located in: output_folder/fireball_name

Most data sources also allow for .kml files to be downloaded. It may be useful to compile all seismic stations into a single kml file for use with Google Earth Pro, and using the circle measure tool to find stations within 150 km of an event.

Sites that stations can be downloaded from:

North Carolina - http://service.ncedc.org/fdsnws/dataselect/1/
South Carolina - https://service.scedc.caltech.edu/fdsnws/dataselect/1/
Turkey - http://tdvm.afad.gov.tr/
Europe - http://udim.koeri.boun.edu.tr/mseed/
Europe - http://eida.gfz-potsdam.de/webdc3/
IRIS - http://service.iris.edu/fdsnws/dataselect/docs/1/builder/

Others may be found by searching for .mseed files online

====================
References
====================
some text


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   docs/About
   docs/Installation
   docs/Atmosphere
   docs/Usage


