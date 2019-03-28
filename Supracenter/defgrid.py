import numpy as np

def defgrid(searchVol, size=64):
    """ Purpose: Spaces the search volume into a grid with a specified number of spacings

    Arguments:
        searchVol: [list] Coordinates to be searched for Supracenters, [x1,x2,y1,y2,z1,z2].
        size: [int] Number of spacings to use in the grid

    Returns:
        x,y,z: [list] Coordinates for the points of the grid that will be searched
        size: [int] The spacing used in the grid
    """

    # Makes spacing based off of given size
    x1 = searchVol[0]
    x2 = searchVol[1]
    y1 = searchVol[2]
    y2 = searchVol[3]
    z1 = searchVol[4]
    z2 = searchVol[5]

    # Fix backwards coordinates
    if x2 < x1:
        x1, x2 = x2, x1

    if y2 < y1:
        y1, y2 = y2, y1
        
    if z2 < z1:
        z1, z2 = z2, z1

    # Division
    dx = (x2 - x1)/(size - 1.0)
    dy = (y2 - y1)/(size - 1.0)
    dz = (z2 - z1)/(size - 1.0)

    x = np.arange(x1, x2 + dx, dx)
    y = np.arange(y1, y2 + dy, dy)
    z = np.arange(z1, z2 + dz, dz)

    return x, y, z, size
