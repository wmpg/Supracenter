# Bolide Acoustic Modelling (BAM)
Computing fireball fragmentation locations and trajectory from seismic or infrasound data.

## Installation

1) This code depends on the WesternMeteorPyLib library. Make sure to install this library first: [https://github.com/dvida/WesternMeteorPyLib](https://github.com/wmpg/WesternMeteorPyLib)

2) In addition to WesternMeteorPyLib's requirements, you will need to install additional libraries. It is possible to install the libraries using pip, but you will have to manaully compile pyhdf, which is not fun. The easiest way to install everything is through conda:

```
conda install -c conda-forge pyhdf simplekml pyopengl pyqtgraph folium
pip install pyqt5
pip install obspy
conda install -c conda-forge netcdf4 pyqtwebengine
pip install pyswarm
```

(Optional) Termcolor is supported for clearer terminal messages, but the program will still work without it:

```
pip install termcolor
```

(Optional) If you need a higher altitude wind model. Temperatures are done through wmpl NRLMSISE. This will require a fortran compiler to run, I use this:
```
conda install -c conda-forge fortran-compiler
pip install hwm93
```


(Optional) If you want to use the signal finder, you'll need this:
```
pip install more-itertools
```


3) Unzip the code into a folder on your computer

## Running a test case

The GUI for the program, which contains all of the features, can be run by:
```
python -m supra.bam
```
In the Anaconda prompt. Make sure you navigate to the folder above supra/ in your Anaconda prompt first.

For further documentation, look for About->Documentation in the GUI.

## Citations

For academic use, please cite the paper:
>[McFadden, L., Brown, P., Vida, D. and Spurný, P., 2021. Fireball characteristics derivable from acoustic data. Journal of Atmospheric and Solar-Terrestrial Physics, 216, p.105587.](https://doi.org/10.1016/j.jastp.2021.105587)
