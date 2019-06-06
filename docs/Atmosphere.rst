.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Atmosphere
**********

=========================
Effect of Weather on Rays
=========================

.. math::
	c_s = \sqrt{\frac{\gamma R T}{M_0}}
.. math::
	c_{eff} = c_s + \hat{n} \cdot \vec{w}

=============================
Atmospheric Data Instructions 
=============================

.. note::
 	Weather data has been tested with recent files, older files may have a different format, in which a custom script to read them is required. See netCDFconv.py for more details. 

Reading the wrong file may cause the program to freeze, such as trying to read UKMO .nc file as an ECMWF .nc file. Killing the program while it is reading the data may corrupt the file (while the message 'Converting weather data. This may take a while...' is in the terminal).  

Constants used in converting weather values can be changed in supra/Fireballs/SeismicTrajectory.ini


===============================
Setting Up Atmospheric Fetching
===============================

Setting up atmospheric fetching is a one-time setup which will allow data to be downloaded from the Copernicus Climate Change Store automatically.

Follow the instructions given here, with the given exceptions below:
https://cds.climate.copernicus.eu/api-how-to

Step 1: Needs to be done locally, since the url and key depend on the users login

Step 2: Proceed as normal

Step 3: Does not need to be done, the Fetch Atmosphere tab does this for you, and will request the data, and display it automatically. 

.. warning::
	Cancelling the program while it is using the weather data may corrupt it! Since these files are so large to download (usually between 1 and 6 GB), it is recommended to keep a copy of the weather file on your computer, just in case!



===================================
UKMO Instructions - Manual Download
===================================
Downloadable files can be found here:
http://data.ceda.ac.uk/badc/ukmo-assim/

An account will need to be created in order to access files

UKMO data is found as a .pp file on the CEDA site

Supracenter only handles .nc files, so the file needs to be converted.

---CONVERT .PP TO .NC-------------------------------------------------
Download xConv1.94 <replace 1.94 with the current version number>
http://cms.ncas.ac.uk/documents/xconv/
and extract xconv1.94 and convsh1.94. Create a symbolic link to convsh: 

Terminal:
ln -s xconv1.94 convsh1.94

In the folder with xconv and convsh, add conv2nc.tcl and runXconv.bat from the xconv_scripts file. Other
example scripts can be found here: http://cms.ncas.ac.uk/documents/xconv/

After this setup, .pp files in the folder can be converted to .nc files with the command:
(Navigate to the folder)

Terminal:
./convsh1.94 ./conv2nc.tcl -i <input file.pp> -o <output file.nc>
----------------------------------------------------------------------

Indicate in the .ini where the .nc file is, and the UKMO data can now be read by Supracenter

Note: killing the program while reading a .nc file WILL corrupt the file

Note: Data conversion may not work for older files, and may need to be read with a script

====================================
ECMWF Instructions - Script Download
====================================
Create an account with ECMWF

Follow the instructions given:
https://software.ecmwf.int/wiki/display/CKB/How+to+download+data+via+the+ECMWF+WebAPI

Follow the instructions given up to before making the script, the script is in mainmenu.py:
https://software.ecmwf.int/wiki/display/CKB/How+to+download+ERA5+data+via+the+ECMWF+Web+API

Change the ECMWF script in mainmenu.py as needed, indicate the hour of the weather profile needed
atm_hour is the hour that will be downloaded, if it is available.
grid_size specifies the spatial resolution of the weather data from ECMWF 

Note: Automatic data retrival may not work for older files

======================================
MERRA-2 Instructions - Script Download
======================================
Create an account: 
https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login

Link the database with your account:
https://disc.gsfc.nasa.gov/earthdata-login

Enter your username and password into fetchMERRA.py

Specify in the .ini where the file is and the name of the file, indicate the hour of the weather profile needed.
atm_hour is the hour that will be downloaded, if it is available.

Note: Automatic data retrival may not work for older files

======================================
MERRA-2 Instructions - Manual Download
======================================
Create an account:
https://wiki.earthdata.nasa.gov/display/EL/How+To+Register+With+Earthdata+Login

Link the database with your accout:
https://disc.gsfc.nasa.gov/earthdata-login

Download the data from the site:
https://disc.gsfc.nasa.gov/datasets/M2I3NVASM_V5.12.4/summary?keywords=%22MERRA-2%22

Go to "subset/get data" and specify what area data is needed. Make sure the file type is HDF
and the variables contain at least "air_temperature", "northward_wind", "eastward_wind", and "mid_layer_heights"

Specify in the .ini where the file is and the name of the file

Note: Data conversion may not work for older files, and may need to be read with a script

==================================
Custom Atmospheric Data
==================================
Create a custom .txt file containing the data, separated by spaces.

There must be at least two rows for the program to be able to run.

The format of the .txt file should be:
--------------------------------------
header
elevation(m) temperature(deg C) wind_speed(m/s) wind_direction(deg from north due east)
--------------------------------------

Specify in the .ini where the file is.

==================================
Isotropic Atmospheric Data
==================================
No file is required for an isotropic profile.

The speed_of_sound variable must be set, and the program will make a custom atmospheric dataset using the given speed of sound.

.. toctree::
   :maxdepth: 2

   





