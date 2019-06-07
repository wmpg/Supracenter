.. Supracenter documentation master file, created by
   sphinx-quickstart on Mon Jun  3 11:53:56 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


INI Files
---------

ini files contain the information for a given fireball event. These are loaded into each part of the program with the parameters of the fireball. For most parts of the program, not all of the variables in the ini files need to be set. For example, in the Supracenter solver, the ballistic search parameters do not need to be filled out.

If the terminal returns an "INI ERROR" while running the code, try setting [General] debug = True, which will print out what variables are being detected by the program.
    
.. toctree::
   :maxdepth: 2

   variables





