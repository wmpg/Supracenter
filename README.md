# Supracenter
Computing fireball fragmentation locations and trajectory from seismic or infrasound data.

## Installation

1) This code depends on the WesternMeteorPyLib library. Make sure to install this library first: [https://github.com/dvida/WesternMeteorPyLib](https://github.com/wmpg/WesternMeteorPyLib)

1) In addition to WesternMeteorPyLib's requirements, you will need to install additional libraries. It is possible to install the libraries using pip, but you will have to manaully compile pyhdf, which is not fun. The easiest way to install everything is through conda:

```
conda install -c conda-forge pyhdf obspy simplekml pyopengl pyqtgraph folium
conda install -c anaconda netcdf4
pip install pyqt5 pyswarm
```

## Running a test case

The GUI for the program, which contains all of the features, can be run by:
```
python -m supra.bam
```

For further documentation, look for About->Documentation in the GUI or if the GUI is not running, look for supra/Fireballs/docs/index.html in the files
