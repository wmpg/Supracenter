import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, poisson, lognorm

FILE_NAME = "F:\\Desktop\\alaska_sens_40.csv"
FILE_NAME2 = "F:\\Desktop\\alaska_sens_constrained.csv"
#FILE_NAME2 = "F:\\Desktop\\alaska_sens_heights.csv"


STATIONS = [[68.640800, -149.572400, "TA-TOLK"], [67.226900, -150.203800, "AK-COLD"]]

R_TOL = 100
YIELD_TOL = 4.184e12
YIELD_MIN = 4.184e6
BINS = 10

def readFile(file_name):
    lats = []
    lons = []
    heights = []
    energies = []
    residuals = []

    with open(file_name, "r+") as f:
        lines = f.readlines()
        for ll, line in enumerate(lines):
            
            if ll == 0:
                continue
            
            line = line.strip("\n")
            data = line.split(",")

            lat = float(data[0])
            lon = float(data[1])
            h = float(data[2])
            e = float(data[4])
            r = float(data[3])

            if not np.isnan(r) and YIELD_MIN < e < YIELD_TOL:
                lats.append(lat)
                lons.append(lon)
                heights.append(h)
                energies.append(e)
                residuals.append(r)


    lats = np.array(lats)
    lons = np.array(lons)
    heights = np.array(heights)
    energies = np.array(energies)
    residuals = np.array(residuals)

    return lats, lons, heights, energies, residuals

lats1, lons1, heights1, energies1, residuals1 = readFile(FILE_NAME)
lats2, lons2, heights2, energies2, residuals2 = readFile(FILE_NAME2)

N1 = len(lats1)
N2 = len(lats2)

plt.ion()
plt.rcParams["figure.figsize"] = (16,9)

r_tol_list = []
mu_list1 = []
err_list1 = []
min_list1 = []
max_list1 = []

mu_list2 = []
err_list2 = []
min_list2 = []
max_list2 = []


while True:
    
    R_TOL -= 2

    keep_ind1 = np.where(residuals1 <= R_TOL)
    rm_ind1 = np.where(residuals1 > R_TOL)

    keep_ind2 = np.where(residuals2 <= R_TOL)
    rm_ind2 = np.where(residuals2 > R_TOL)

    lats1_filt = lats1[keep_ind1]
    lons1_filt = lons1[keep_ind1]
    heights1_filt = heights1[keep_ind1]
    energies1_filt = energies1[keep_ind1]/4.184e6
    residuals1_filt = residuals1[keep_ind1]

    lats2_filt = lats2[keep_ind2]
    lons2_filt = lons2[keep_ind2]
    heights2_filt = heights2[keep_ind2]
    energies2_filt = energies2[keep_ind2]/4.184e6
    residuals2_filt = residuals2[keep_ind2]

    n1 = len(residuals1_filt)
    n2 = len(residuals2_filt)        

    if n1 < 10 or n2 < 10:
        break
    

    plt.clf()

    plt.subplot(231)


    plt.scatter(np.log10(np.array(energies1[rm_ind1])/4.186e6), residuals1[rm_ind1], c='r')
    plt.scatter(np.log10(np.array(energies1_filt)), residuals1_filt, c='k', label="File 1", alpha=0.3)

    plt.scatter(np.log10(np.array(energies2[rm_ind2])/4.186e6), residuals2[rm_ind2], c='r')
    plt.scatter(np.log10(np.array(energies2_filt)), residuals2_filt, c='g', label="File 2", alpha=0.3)

    plt.axhline(R_TOL)
    plt.semilogx()
    plt.xlabel("log10(Energy) [kT TNT]")
    plt.ylabel("Error in Solution [s]")

    plt.legend()

    plt.xlim([np.log10(min(energies1)/4.186e6), np.log10(max(energies1)/4.186e6)])
    plt.ylim([0, max(residuals1)])

    plt.subplot(232)

    plt.tricontour(lons1_filt, lats1_filt, residuals1_filt, linewidths=0.5, colors='w')
    cntr = plt.tricontourf(lons1_filt, lats1_filt, residuals1_filt, cmap="viridis")
    plt.colorbar(cntr)

    for stat in STATIONS:
        plt.scatter(stat[1], stat[0], marker="^", c="k")
        plt.text(stat[1], stat[0], stat[2])


    plt.xlabel("Longitude [E]")
    plt.ylabel("Latitude [N]")
    # plt.xlim([-149.5, -148.5])
    # plt.ylim([67.5, 68.5])

    plt.subplot(212)

    plt.hist(np.log10(energies1_filt), density=True, bins=int(n1/BINS), alpha=0.5, color='b', label="File 1")
    plt.hist(np.log10(energies2_filt), density=True, bins=int(n2/BINS), alpha=0.5, color='g', label="File 2")


    plt.xlim([np.log10(min(energies1)/4.186e6), np.log10(max(energies1)/4.186e6)])
    # plt.ylim([0, 1])

    def getHist(energy):
        mu, std = norm.fit(np.log10(energy))
        xmin, xmax = plt.xlim()
        x = np.linspace(xmin, xmax, 100)
        p = norm.pdf(x, mu, std)

        return x, p, mu, std

    x1, p1, mu1, std1 = getHist(energies1_filt)
    x2, p2, mu2, std2 = getHist(energies2_filt)  

    # plt.plot(x, lognorm.pdf(x, mu), label='Log Norm')
    plt.plot(x1, p1, 'b', linewidth=2, label="File 1: Normal Distribution: Mean {:.2f} ({:.2f}) kT TNT".format(10**mu1, 10**std1))
    plt.plot(x2, p2, 'g', linewidth=2, label="File 2: Normal Distribution: Mean {:.2f} ({:.2f}) kT TNT".format(10**mu2, 10**std2))
    plt.title("N1 = {:}, N2 = {:}, R_TOL = {:}".format(n1, n2, R_TOL))
    plt.ylabel("Normalized Number of Observations")
    plt.xlabel("log10 Energy [kT TNT]")
    plt.legend()


    plt.subplot(233)

    SIGMA = 3

    r_tol_list.append(R_TOL)
    mu_list1.append(10**mu1)
    err_list1.append(SIGMA*10**std1)
    min_list1.append(max([min(energies1_filt), 10**mu1 - 3*10**std1]))
    max_list1.append(min([max(energies1_filt), 10**mu1 + 3*10**std1]))

    mu_list2.append(10**mu1)
    err_list2.append(SIGMA*10**std2)
    min_list2.append(max([min(energies2_filt), 10**mu2 - 3*10**std2]))
    max_list2.append(min([max(energies2_filt), 10**mu2 + 3*10**std2]))

    plt.scatter(r_tol_list, mu_list1, c='k', alpha=0.3)
    plt.errorbar(r_tol_list, mu_list1, yerr=err_list1, fmt="o", c='k', capsize=5, label="File 1: 3 Std Dev")

    plt.scatter(r_tol_list, mu_list2, c='g', alpha=0.3)
    plt.errorbar(r_tol_list, mu_list2, yerr=err_list2, fmt="o", c='g', capsize=5, label="File 2: 3 Std Dev")

    # plt.errorbar(r_tol_list, mu_list, yerr=[min_list, max_list], fmt="o", c='r', capsize=5)
    plt.xlabel("Allowed Tolerance [s]")
    plt.ylabel("Energy Estimate [kT TNT]")
    plt.gca().invert_xaxis()
    plt.legend()

    plt.draw()
    plt.pause(0.0001)
    
plt.waitforbuttonpress()