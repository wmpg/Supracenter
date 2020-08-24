import numpy as np

from supra.Utils.Classes import Position
from supra.Supracenter.anglescanrev import anglescanrev


def propegateBackwards(ref_pos, stn, bam, offset=0):
    
    S = stn.metadata.position
    
    S.pos_loc(ref_pos)

    # Initial guess for sounding
    sounding, _ = bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, 50000])

    D = []
    offset = 0
    for pol in range(len(stn.polarization.azimuth)):
        min_az = stn.polarization.azimuth[pol] - stn.polarization.azimuth_error[pol]
        max_az = stn.polarization.azimuth[pol] + stn.polarization.azimuth_error[pol]

        for azimuth in np.linspace(min_az, max_az, 10):
            for zenith in np.linspace(1, 89, 50):
                
                # T - expected final arrival, with bad sounding
                # Recalculate winds
                # D - real, expected final arrival
                # This is overkill, atmosphere won't change that much

                T_pos = anglescanrev(S.xyz, (azimuth + offset)%360, zenith, sounding, wind=True)
                T = Position(0, 0, 0)
                T.x = T_pos[0]
                T.y = T_pos[1]
                T.z = T_pos[2]
                T.pos_geo(ref_pos)

                try:
                    sounding_plus, perts = bam.atmos.getSounding(lat=[S.lat, T.lat], lon=[S.lon, T.lon], heights=[S.elev, 50000])
                except ValueError:
                    sounding_plus, perts = bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, 50000])
                
                nom_data = [None]*len(perts+1)
                nom_data[0] = anglescanrev(S.xyz, (azimuth + offset)%360, zenith, sounding_plus, wind=True, trace=True)
                for pp, pert in enumerate(perts):
                    nom_data[pp] = anglescanrev(S.xyz, (azimuth + offset)%360, zenith, pert, wind=True, trace=True)


                # Repeat 180 deg away

                T_pos = anglescanrev(S.xyz, (azimuth + offset + 180)%360, zenith, sounding, wind=True)
                T = Position(0, 0, 0)
                T.x = T_pos[0]
                T.y = T_pos[1]
                T.z = T_pos[2]
                T.pos_geo(ref_pos)

                try:
                    sounding_plus, perts = bam.atmos.getSounding(lat=[S.lat, T.lat], lon=[S.lon, T.lon], heights=[S.elev, 50000])
                except ValueError:
                    sounding_plus, perts = bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, 50000])
                
                nom_data_rev = [None]*len(perts+1)
                nom_data_rev[0] = anglescanrev(S.xyz, (azimuth + offset + 180)%360, zenith, sounding_plus, wind=True, trace=True)
                for pp, pert in enumerate(perts):
                    nom_data_rev[pp] = anglescanrev(S.xyz, (azimuth + offset + 180)%360, zenith, pert, wind=True, trace=True)



                for ii in range(len(nom_data)):
                    for line in nom_data[ii]:
                        line[3] -= stn.polarization.time[pol]
                        D.append(line)

                    for line in nom_data_rev[ii]:
                        line[3] -= stn.polarization.time[pol]
                        D.append(line)

    return D

