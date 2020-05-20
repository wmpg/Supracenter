.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


All Times Files
---------------

These files contain the estimated arrival times for each perturbation at each station for all fragmentation and ballistic arrivals. This .npy file is used to quickly save and load data and prevent long loading times. The data can be called from a numpy array with the data loaded in as follows:

[perturbation #, station #, 0 - ballistic/ 1 - fragmentation, frag number (0 if using ballistic)]

The intended use of this is to save estimates for arrival times that take a long time to load, such as perturbation calculations. It is recommended that each .npy file which contains perturbations is saved in a proper place so that it may be reused.

All Times files are automatically created by the Make Picks window when the stations are loaded in, under **all_pick_times.npy** in your current working directory, but a file may be exported by the Make Picks window by selecting **Export All Times File**.

All Times files may be loaded in through the ini file being used.

Theoretical All Times Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^

An all times file created by the TheoreticalPicks code has two extra entries (zenith and velocity), and is missing the frag number tag. Therefore, it is not compatable with the rest of the program.

[perturbation, station, 0 - ballistic/ 1 - fragmentation, zeniths, velocity]


.. toctree::
   :maxdepth: 2

   





