import numpy as np

from supra.Atmosphere.Pressure import *
from supra.Yields.YieldFuncs import *
from supra.Yields.YieldCalcs import *



# Problem 7-1
print("EXAMPLE 7-1")
W_0 = 4.184e6
W = 100*W_0

W_W_0 = W/W_0

H = 10000
R = 15

# pressure at H
P = estPressure(H)

# Pressure at sea level
P_0 = estPressure(0)

print("Pressure at {:.2f} km: {:.2f} Pa ({:.2f} mb)".format(H/1000, P, pascal2Millibar(P)))

# Transmission Factor
T = transmissionFactor(H, mode="alt")

print("Transmission Factor: {:.3f}".format(T))

# Scaled distance - note, P_P_0 is 1 since the explosion happens 
# at the same height as the reciever
R_0 = scaledDistance(R, W_W_0, P/P, f_T=T)

print("Scaled Distance: {:.2f} m".format(R_0))

# Overpressure ratio of the scaled distance
p_rat = chemOverpressureRatio(R_0)

print("Peak Overpressure Ratio at Sea Level: {:.2f}".format(p_rat))
print("Peak Overpressure at Height: {:.2f} Pa ({:.2f} mb)".format(p_rat*P, pascal2Millibar(p_rat*P)))

print("")
print("EXAMPLE 7-2")
# Problem 7-2
W = Yield(3000*4.184e12)

H1 = 10000
P1 = estPressure(H1)

H2 = 4000
P2 = estPressure(H2)

R = 10000

P1 = 26500
P2 = 61700

print(W)
print("Height 1 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H1/1000, P1, pascal2Millibar(P1)))
print("Height 2 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H2/1000, P2, pascal2Millibar(P2)))

T1 = transmissionFactor(H1, mode="mean")
T2 = transmissionFactor(H2, mode="mean")
T = (H1*T1 - H2*T2)/(H1 - H2)
T = 0.784

print("Transmission Factor: H1 = {:.3f}, H2 = {:.3f}, Total = {:.3f}".format(T1, T2, T))

Z = T*R/(W.nuclear)**(1/3)

print("Scaled Distance: {:.2f} m".format(Z))
print("Actual Distance: {:.2f} m".format(R))

# Convert to overpressure
# Overpressure ratio of the scaled distance (this needs to be nuc to work for this question, which isn't complete yet)
p_rat = chemOverpressureRatio(Z)
p_rat = 0.230
print("Overpressure Ratio {:.3f}".format(p_rat))
print("Overpressure: {:.2f} Pa ({:.2f} mb)".format(p_rat*P2, pascal2Millibar(p_rat*P2)))

print("")
print("EXAMPLE 7-3")

R = 3.75
del_p = millibar2Pascal(1225)
P = estPressure(0)

p_rat = del_p/P
print("Overpressure Ratio: {:.3f}".format(p_rat))

Z = chemOverpressureRatioInv(p_rat)
print("Actual Distance: {:.2f} m".format(R))
print("Scaled Distance: {:.2f} m".format(Z))

W = overpressure2YieldKG(del_p, 0, 0, R)

print("")
print("EXAMPLE 7-2+")
# Problem 7-2
W = Yield(3000*4.184e6)

H1 = 10000
P1 = estPressure(H1)

H2 = 4000
P2 = estPressure(H2)

R = 10000

P1 = 26500
P2 = 61700

yield2OverpressureKG(W.j, H1, H2, R)
print("")