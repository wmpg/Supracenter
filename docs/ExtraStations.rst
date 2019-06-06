.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Extra Stations
**************

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


.. toctree::
   :maxdepth: 2
   :caption: Contents:





