# Supracenter
Computing fireball fragmentation locations and trajectory from seismic or infrasound data.

## Installation

1) This code depends on the WesternMeteorPyLib library. Make sure to install this library first: [https://github.com/dvida/WesternMeteorPyLib](https://github.com/dvida/WesternMeteorPyLib])

1) In addition to WesternMeteorPyLib's requirements, you will need to install additional libraries. It is possible to install the libraries using pip, but you will have to manaully compile pyhdf, which is not fun. The easiest way to install everything is through conda:

```
conda install -c conda-forge pyhdf obspy simplekml
conda install -c anaconda netcdf4
pip install pyqt5 pyswarm
```

## Running a test case

Create an .ini file using either the example .ini located in ```supra/Fireballs/docs/example_ini.ini``` or by creating a custom .ini file using the information given in the user manual in ```supra/Fireballs/docs/User Manual.txt```

First, the waveform data must be obtained by running:
```
python -m wmpl.Fireballs.GetIRISData <location of .ini file>
```
with the variable get_data = True

To analyze the waveform data and make arrival time picks:
```
python -m wmpl.Fireballs.MakeIRISPicks <location of .ini file>
```

If picks have been created, the Supracenter module can be run by:
```
python -m wmpl.Supracenter.mainmenu <location of .ini file>
```

and the Trajectory solution can be run by:
```
python -m wmpl.Fireballs.SeismicTrajectory <location of .ini file>
```

## More documentation

You can find more documentation under ```supra/Fireballs/docs``` and ```supra/Fireballs/Supracenter```.
