
from supra.Atmosphere.Pressure import estPressure

from supra.Yields.YieldFuncs import *


#B1.1 10 kT NE at 3 km altitude
print("########## B1.1: 10 kT NE at 3 km altitude")

del_p_0 = 1
R_0 = 1

W = Yield(10*4.184e12)

H_1 = 3000 #m
H_0 = 0

print(W)

# pressure at ground
P_0 = estPressure(H_0)

# pressure at 3 km
P = estPressure(H_1)

print("Pressure at ground: {:.2f} kPa".format(P_0/1000))
print("Pressure at 3 km  : {:.2f} kPa".format(P/1000))

P_P_0 = getPP0(H_1)

# scaled overpressure
del_p_1 = scaledOverpressure(del_p_0, P_P_0)

# scaled distance
R_1 = scaledDistance(R_0, W.nuclear, P_P_0)
# R_1 = (W_1/W_0)**(1/3)*(P_0/P)**(1/3)*R_0

print("Scaled Pressure: {:.2f} Pa".format(del_p_1))
print("Scaled Distance: {:.2f} m".format(R_1))

# Scaled Height of Burst (HOB)
HOB = H_1/W.ysf

print("Scaled HOB: {:.2f} m".format(HOB))

print("")

print("########## B1.2: 1000 kg HE Surface Burst at Sea Level")

# 4000 kg HE free air burst (shouldn't be NE like in the text)
W = Yield(4.4e-3*4.184e12)

print(W)

print("Yield Scaling Factor: {:.3f}".format(W.ysf))


print("")

#B1.3 100 kg HE at 20 m HOB
print("########## B1.3: 100 kg HE at 20 m HOB")


H_1 = 20 #m
H_0 = 0 

W = Yield(100*4.184e6)

# scaled HOB
HOB = scaledHOB(H_1, W.chem)

print("Scaled HOB: {:.2f} m".format(HOB))

# Scaled HOB factor
HOB_scale = r02HOB(HOB)

print("Yield Scaling Factor: {:.2f}".format(HOB_scale))

# Equivalent airburst yield
air_W = HOB_scale*W.chem

print("Equivalent Airburst Yield: {:.2f} kg HE".format(air_W))

# convert from kg HE to kT NE
air_W_NE = kgHE2kTNE(air_W)

print("Equivalent Airburst Yield: {:.2e} kT NE".format(air_W_NE))

# yield scaling factor
ysf = yieldScalingFactor(air_W_NE)

print("Yield Scaling Factor: {:.3f}".format(ysf))

R = 1000

# scaled distance
R_0 = scaledDistance(R, air_W_NE, 1)

print("Scaled Distance: {:.2f} m".format(R_0))

# find overpressure
del_p = overpressureDistance(R, air_W_NE, 1, mode="kTNE")

print("Overpressure: {:.2f} Pa".format(del_p))



