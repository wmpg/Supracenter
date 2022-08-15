import numpy as np


from supra.Atmosphere.Pressure import *
from supra.Yields.YieldFuncs import *
from supra.Yields.YieldCalcs import *


print("Section 12.3 Examples of Scaling")



W_0 = 1 # in pounds
W = 1000 # in pounds

d = 3.61
print("We know that {:} ft away from a {:} pound charge produces 60 psi".format(d, W_0))
print("Therefore, a {:.2f} pound charge will produce the same overpressure at {:.2f} ft".format(W, (W/W_0)**(1/3)*d))


print("To find the range at which 110 kPa occurs for a 1 g charge")
f = (1/454)**(1/3)
print("Our radius is multiplied by a factor of {:.3f}".format(f))
# 1 pound = 454 grams
print("So if 110 kPa occurs at 1.91 m for a 1 pound charge, it will occur at {:.4f} m for a 1 g charge".format(f*1.91))


print("")
print("Atmospheric Scaling")

H1 = 6500/3.281 # meters
H2 = 0

P1 = estPressure(H1)
P2 = estPressure(H2)

print("Blast from a height of {:.2f} km ({:.2f} Pa) to the ground ({:.2f} Pa)".format(H1/1000, P1, P2))
print("Atmospheric Pressure Ratio: {:.3f}".format(P1/P2))

print("The distance to 15.4 psi is 6.4 ft from the table (interpolated)")
print("This is also the distance to 12 psi at 6500 ft")


