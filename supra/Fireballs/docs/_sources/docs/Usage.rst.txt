.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. # with overline, for parts
.. * with overline, for chapters
.. =, for sections
.. -, for subsections
.. ^, for subsubsections
.. ", for paragraphs

#####
Usage
#####

Main Program
============

To open the program, navigate to your Supracenter folder and type the following command::

    python -m supra.Fireballs.SolutionGUI

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   inibuilder
   stations
   makepicks
   fetchatmosphere
   atmosphereviewer
   pickviewer
   supracentersearch
   supracentermanual
   seismictraj

Controls:
* Fullscreen Toggle: F11
* INI Toolbar Toggle: V
* Save: CTRL+S
* Save As: CTRL+SHIFT+S
* Quit: CTRL+Q

Other Functions
===============

These functions may be used independently of the program, but do not have an interface.

TrajectoryInterp.py
-------------------
Will find *N* points along a trajectory from either a start and end point, or a start point and the azimuth and zenith angle.

TheoreticalPicks.py
-------------------
Used in combination with allTimesRead.py, can model theoretical ballistic arrivals of meteors with different velocity and zenith angles. 






