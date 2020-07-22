from supra.Utils.Classes import Position

import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.cyscan2angle import cyscan

def finalanglecheck(bam, traj, D, angle):

    points = traj.trajInterp2(div=100)

    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)

    azimuths = []

    for pt in points:

        S = Position(pt[0], pt[1], pt[2])
        
        S.pos_loc(ref_pos)
        D.pos_loc(ref_pos)
        
        sounding, _ = bam.atmos.getSounding(lat=[S.lat, D.lat], lon=[S.lon, D.lon], heights=[S.elev, D.elev])
        
        az = cyscan(S.xyz, D.xyz, sounding)

        azimuths.append(az%180)


    best_az_indx = np.nanargmin(np.abs(azimuths - angle))

    err = np.abs(azimuths[best_az_indx] - angle)

    return points[best_az_indx], err

