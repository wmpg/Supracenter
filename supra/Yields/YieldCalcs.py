
import numpy as np

from supra.Atmosphere.Pressure import *
from supra.Yields.YieldFuncs import *
from supra.Utils.pso import pso


def getWW0(W, mode="kgHE"):

    if mode == "kgHE":
        W_0 = 4.184e6

    elif mode == "kTNE":
        W_0 = 4.184e12

    else:
        print("Unknown mode: {:}".format(mode))
        return None

    W_W_0 = W/W_0

    return W_W_0



def yield2Overpressure(W, R, HOB, mode="kgHE", verbose=False):

    # pressure at ground
    P_0 = estPressure(0)

    # pressure at HOB
    P = estPressure(HOB)

    W_W_0 = getWW0(W, mode=mode)

    P_P_0 = getPP0(HOB)

    # scaled HOB
    HOB_s = scaledHOB(HOB, W_W_0)

    # Scaled HOB factor
    HOB_factor = r02HOB(HOB_s)

    # Equivalent airburst yield
    air_W = HOB_factor*W_W_0

    if mode == "kgHE":
        # convert from kg HE to kT NE
        air_W_NE = kgHE2kTNE(air_W)

    # find overpressure
    del_p = overpressureDistance(R, air_W_NE, P_P_0, mode="kTNE")

    # scaled distance
    R_0 = scaledDistance(R, air_W_NE, P_P_0)
        
    if verbose:
        print("####### RESULTS #######")
        print("Estimation of {:.2E} J burst at {:.2f} km".format(W, HOB/1000))
        print("Range to Source {:.2f} km".format(R/1000))
        print("Pressure at ground: {:.2f} kPa".format(P_0/1000))
        print("Pressure at 3 km  : {:.2f} kPa".format(P/1000))
        print("Pressure Ratio:     {:.2f}".format(P_P_0))
        print("Yield conversion {:.2f} {:}".format(W_W_0, mode))
        print("Scaled HOB: {:.2f} m".format(HOB_s))
        print("HOB Scaling Factor: {:.2f}".format(HOB_factor))
        print("Equivalent Airburst Yield: {:.2f} kg HE".format(air_W))

        if mode == "kgHE":
            print("Equivalent Airburst Yield: {:.2e} kT NE".format(air_W_NE))

        print("Overpressure: {:.2f} Pa".format(del_p))
        print("Scaled Distance: {:.2f} m".format(R_0))

    return del_p

def yield2OverpressureKG(W, H1, H2, R):

    ### Only working work chem
    W = Yield(W)

    P1 = estPressure(H1)
    P2 = estPressure(H2)


    print(W)
    print("Height 1 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H1/1000, P1, pascal2Millibar(P1)))
    print("Height 2 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H2/1000, P2, pascal2Millibar(P2)))

    T = 1

    Z = T*R/(W.chem)**(1/3)

    print("Scaled Distance: {:.2f} m".format(Z))
    print("Actual Distance: {:.2f} m".format(R))

    # Convert to overpressure
    # Overpressure ratio of the scaled distance
    p_rat = chemOverpressureRatio(Z)


    print("Overpressure Ratio {:.3f}".format(p_rat))
    print("Overpressure: {:.2f} Pa ({:.2f} mb)".format(p_rat*P2, pascal2Millibar(p_rat*P2)))

    return p_rat*P2

def yield2OverpressureKGNuc(W, H1, H2, R):

    W = Yield(W)

    P1 = estPressure(H1)
    P2 = estPressure(H2)


    print(W)
    print("Height 1 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H1/1000, P1, pascal2Millibar(P1)))
    print("Height 2 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H2/1000, P2, pascal2Millibar(P2)))

    T = 1

    Z = T*R/(W.nuclear)**(1/3)

    print("Scaled Distance: {:.2f} m".format(Z))
    print("Actual Distance: {:.2f} m".format(R))

    p_rat = nucOverpressureFromScaledDistance(Z)

    print("Overpressure Ratio {:.3f}".format(p_rat))
    print("Overpressure: {:.2f} Pa ({:.2f} mb)".format(p_rat*P2, pascal2Millibar(p_rat*P2)))

    return p_rat*P2

def overpressure2YieldKGNuc(del_p, H1, H2, R):

    P1 = estPressure(H1)
    P2 = estPressure(H2)

    print("Height 1 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H1/1000, P1, pascal2Millibar(P1)))
    print("Height 2 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H2/1000, P2, pascal2Millibar(P2)))

    T = 1

    p_rat = del_p/P2

    print("Overpressure Ratio: {:.3f}".format(p_rat))

    Z = nucScaledDistanceFromOverpressure(p_rat)

    print("Scaled Distance: {:.2e} m".format(Z))
    print("Actual Distnace: {:.2e} m".format(R))

    W = Yield((T*R/Z)**3*4.186e12)

    print(W)

    return W.j


def overpressure2YieldKG(del_p, H1, H2, R):

    P1 = estPressure(H1)
    P2 = estPressure(H2)

    print("Height 1 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H1/1000, P1, pascal2Millibar(P1)))
    print("Height 2 {:.2f} km ({:.2f} Pa = {:.2f} mb)".format(H2/1000, P2, pascal2Millibar(P2)))

    T = transmissionFactor(H1)
    print("Transmission Factor: {:.3f}".format(T))

    
    print("Observed Overpressure: {:.2f} Pa".format(del_p))

    p_rat = del_p/P2
    print("Overpressure Ratio: {:.2e}".format(p_rat))

    Z = chemOverpressureRatioInv(p_rat)

    print("Scaled Distance: {:.2e} m".format(Z))
    print("Actual Distnace: {:.2e} m".format(R))

    W = Yield((R/Z/T)**3*4.186e6)

    print(W)

    return W.chem

def overpressure2Yield(del_p, R, H, verbose=True):

    del_p_orig = del_p

    # pressure at ground
    P_0 = estPressure(0)

    # pressure at H
    P = estPressure(H)

    P_P_0 = getPP0(H)

    # # scale overpressure
    # del_p = scaledOverpressure(del_p_orig, P_P_0)

    # find unscaled yield
    air_W_NE = distance2Overpressure(R, del_p, P_P_0, mode="kgHE")

    W = Yield(air_W_NE*4.184e6)

    # Transmission Factor
    T = transmissionFactor(H, mode="mean")

    # scaled distance
    Z = scaledDistance(R, W.chem, P_P_0, f_T=T)

    # def searchLoop(W_W_0, *loopargs):

    #     HOB, air_W = loopargs

    #     W_W_0 = 10**W_W_0

    #     # scaled HOB
    #     HOB_s = scaledHOB(H, W_W_0)

    #     # Scaled HOB factor
    #     HOB_factor = r02HOB(HOB_s)

    #     # Equivalent airburst yield
    #     air_W_test = HOB_factor*W_W_0

    #     error = np.abs(air_W_test - air_W)

    #     return error


    # f_opt, x_opt = pso(searchLoop, [-3], [8], \
    #     args=[HOB, air_W])

    # W_W_0 = 10**f_opt[0]

    # # scaled HOB
    # HOB_s = scaledHOB(HOB, W_W_0)

    # # Scaled HOB factor
    # HOB_factor = r02HOB(HOB_s)

    if verbose:
        print("####### RESULTS #######")
        print("Estimation of {:.2f} Pa overpressure at {:.2f} km".format(del_p_orig, H/1000))
        print("Scaled Overpressure: {:.2f} Pa".format(del_p))
        print("Scaled Distance: {:.2f} km".format(Z/1000))
        print("Range to Source {:.2f} km".format(R/1000))
        print("Pressure at ground: {:.2f} kPa".format(P_0/1000))
        print("Pressure at {:.1f} km  : {:.2f} kPa".format(H/1000, P/1000))
        print("Pressure Ratio:     {:.2f}".format(P_P_0))
        print("Transmission Factor:{:.3f}".format(T))
        print(W)
        # print("Scaled HOB: {:.2f} m".format(HOB_s))
        # print("HOB Scaling Factor: {:.2f}".format(HOB_factor))


    return W.j



if __name__ == "__main__":

    # del_p = millibar2Pascal(1225)
    H1 = 0
    H2 = 40000
    R = 56569

    op = np.logspace(-2, 2)
    p_rat = op

    W_list = []
    for p in p_rat:
        W = overpressure2YieldKG(p, H1, H2, R)
        W_list.append(W)


    W_list = np.array(W_list)


    plt.plot(p_rat, W_list)
    plt.semilogy()
    plt.xlabel("Overpressure [Pa]")
    plt.ylabel("Explosive Yield [kg TNT HE]")
    plt.show()

    # p_rat = 1
    # R = np.logspace(1, 5)

    # W_list = []
    # for r in R:
    #     W = overpressure2YieldKG(p_rat, H1, H2, r)
    #     W_list.append(W)

    # W_list = np.array(W_list)


    # plt.plot(R, W_list/4.184e12)
    # plt.loglog()
    # plt.xlabel("Range")
    # plt.ylabel("Explosive Yield [kT TNT]")
    # plt.show()


    # #Z = np.array([10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 270, 280, 290, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 800, 850, 900, 950, 1000, 1500, 2000, 3000, 4000, 5000])
    # # p_rat = np.array([3200, 970, 420, 220, 128, 82.9, 57, 41.2, 31, 24, 19.1, 15.5, 12.9, 10.8, 9.22, 6.93, 5.41, 4.35, 3.58, 3, 2.562, 2.215, 1.937, 1.711, 1.524, 1.369, 1.237, 1.125, 1.075, 1.028, 0.985, 0.945, 0.907, 0.871, 0.838, 0.807, 0.778, 0.75, 0.7, 0.655, 0.614, 0.557, 0.5, 0.439, 0.389, 0.348, 0.314, 0.285, 0.261, 0.239, 0.221, 0.205, 0.191, 0.178, 0.167, 0.157, 0.148, 0.140, 0.133, 0.126, 0.114, 0.104, 0.096, 0.088, 0.082, 0.046, 0.032, 0.019, 0.014, 0.011])
    # Z = np.logspace(-1, 6, 100)
    # # r = 1/Z

    # # z = np.polyfit(r, p_rat, 24)
    # # func = np.poly1d(z)

    # # np.save("supra\\Yields\\examples\\nuc_scaled_curve_fit_inv.npy", func)
    
    # # plt.scatter(Z, p_rat)

    # p_rat = []
    # z_fit = []
    # p_fit = []
    # for zz in Z:
    #     p_rat.append(chemOverpressureRatio(zz))
    #     if zz >= 10 and zz <= 200:
    #         z_fit.append(zz)
    #         p_fit.append(chemOverpressureRatio(zz))

    # z_fit = np.array(z_fit)
    # p_fit = np.array(p_fit)

    # plt.scatter(Z, np.array(p_rat))
    # plt.scatter(z_fit, p_fit)

    # r = np.log10(z_fit)
    # s = np.log10(p_fit)
    # z = np.polyfit(r, s, 1)
    # func = np.poly1d(z)

    # plt.plot(Z, 10**func(np.log10(Z)))

    # def chem_func(Z):

    #     p_Pa = 808*(1 + (Z/4.5)**2)/(1 + (Z/0.048)**2)**0.5/(1 + (Z/0.32)**2)**0.5/(1 + (Z/1.35)**2)**0.5
    #     return p_Pa

    # plt.plot(Z, chem_func(Z))

    # np.save("supra\\Yields\\nuc_scaled_curve_fit_extended.npy", func)

    # plt.scatter([100, 125, 150, 200, 250, 300, 400, 500], [0.008, 0.007, 0.005, 0.004, 0.003, 0.003, 0.002, 0.002])
    # plt.ylabel("Pressure Ratio")
    # plt.xlabel("Scaled Distance [m from a 1kT NE explosion]")
    # plt.loglog()
    # plt.show()
    # # W = 1e6*4.184e6
    # HOB = 40000
    # R = 100000
    # del_p = 1.17


    # W_W_0 = 1e5
    # P_P_0 = getPP0(HOB)

    # T = transmissionFactor(HOB)

    # # del_p = yield2Overpressure(W_W_0*4.184e6, R, HOB)
    # W_op = overpressure2Yield(del_p, R, HOB)
    # W_W_0 = W_op / 4.184e6
    # Z = scaledDistance(R, W_W_0, P_P_0, f_T=T)
    # W = yieldFromScaledDistance(Z, R, P_P_0, f_T=T, mode="kgHE")

    # print("Source Height   = {:.2f} km".format(HOB/1000))
    # print("Actual Distance = {:.2f} km".format(R/1000))
    # print("Scaled Distance = {:.2f} m".format(Z))
    # print("Expected Overpressure = {:.2f} Pa".format(del_p))
    # print("Transmission Factor = {:.3f}".format(T))
    # print("Yield           = {:.2f} kg HE".format(W/4.184e6))
    # print("Del_p Yield     = {:.2f} kg HE".format(W_op/4.184e6))
    # print("Initial Yield   = {:.2f} kg HE".format(W_W_0))

    # yield2Overpressure(W, R, HOB)
    # print("")
    