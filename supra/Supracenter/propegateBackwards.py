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
    for zenith in np.linspace(91, 179, 25):
        
        # T - expected final arrival, with bad sounding
        # Recalculate winds
        # D - real, expected final arrival

        T_pos = anglescanrev(S.xyz, (stn.polarization.azimuth + offset)%360, zenith, sounding, wind=True)
        T = Position(0, 0, 0)
        T.x = T_pos[0]
        T.y = T_pos[1]
        T.z = T_pos[2]
        T.pos_geo(ref_pos)

        try:
            sounding_plus, _ = bam.atmos.getSounding(lat=[S.lat, T.lat], lon=[S.lon, T.lon], heights=[S.elev, 50000])
        except ValueError:
            sounding_plus, _ = bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, 50000])
        
        nom_data = anglescanrev(S.xyz, (stn.polarization.azimuth + offset)%360, zenith, sounding_plus, wind=True, trace=True)

        # Repeat 180 deg away

        T_pos = anglescanrev(S.xyz, (stn.polarization.azimuth + offset + 180)%360, zenith, sounding, wind=True)
        T = Position(0, 0, 0)
        T.x = T_pos[0]
        T.y = T_pos[1]
        T.z = T_pos[2]
        T.pos_geo(ref_pos)

        try:
            sounding_plus, _ = bam.atmos.getSounding(lat=[S.lat, T.lat], lon=[S.lon, T.lon], heights=[S.elev, 50000])
        except ValueError:
            sounding_plus, _ = bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, 50000])
        
        nom_data_rev = anglescanrev(S.xyz, (stn.polarization.azimuth + offset + 180)%360, zenith, sounding_plus, wind=True, trace=True)

        for line in nom_data:
            D.append(line)

        for line in nom_data_rev:
            D.append(line)

    return D

