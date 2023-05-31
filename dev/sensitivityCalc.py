import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm, poisson, lognorm
from scipy.optimize import curve_fit

#FILE_NAME = "F:\\Desktop\\alaska_sens_40.csv"
#FILE_NAME = "F:\\Desktop\\alaska_sens_constrained.csv"
#FILE_NAME = "F:\\Desktop\\alaska_sens_heights.csv"
#FILE_NAME = "F:\\Desktop\\alaska_ANSI_KG85_winds.csv"
# FILE_NAME = "F:\\Desktop\\alaska_07_19_22_h.csv"
#FILE_NAME = "F:\\Desktop\\romania_sens_winds_10_grid.csv"
FILE_NAME = "F:\\Desktop\\NY_03_07_2023.csv"
# FILE_NAME = "F:\\Desktop\\NZ_07_28_22.csv"
# FILE_NAME = "F:\\Desktop\\Nord_08_02_22.csv"


CONTOUR = True
BIMODAL = False

#STATIONS = [[68.640800, -149.572400, "TA-TOLK"], [67.226900, -150.203800, "AK-COLD"]]
STATIONS = []
R_TOL = R_TOL_INIT = 100
YIELD_TOL = 4.184e14
YIELD_MIN = 4.184e6
BINS = 10

STATS_TOL = 2


def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1) + gauss(x,mu2,sigma2,A2)

def readFile(file_name):
    print("Reading File...")
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
            stats = float(data[5])

            if not np.isnan(r) and YIELD_MIN < e < YIELD_TOL and stats >= STATS_TOL:
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

def findBestSol(lats1, lons1, heights1, energies1, residuals1):

    best_res = np.nanargmin(residuals1)

    return [lats1[best_res], lons1[best_res], heights1[best_res], energies1[best_res], residuals1[best_res]]

def findAvgSol(lat, lon, h, e, res):

    return [np.mean(lat), np.mean(lon), np.mean(h), np.mean(e), np.mean(res)]


lats1, lons1, heights1, energies1, residuals1 = readFile(FILE_NAME)
print("Initializing...")
N = len(lats1)

plt.ion()
plt.rcParams["figure.figsize"] = (16,9)

r_tol_list = []
mu_list = []
err_list = []
min_list = []
max_list = []

best_sol = findBestSol(lats1, lons1, heights1, energies1, residuals1)

print("Restricting to {:} stations".format(STATS_TOL))

print("Starting Loop")
while True:
    
    R_TOL -= 2

    print("Current R_TOL: {:}".format(R_TOL))

    keep_ind = np.where(residuals1 <= R_TOL)
    rm_ind = np.where(residuals1 > R_TOL)


    lats1_filt = lats1[keep_ind]
    lons1_filt = lons1[keep_ind]
    heights1_filt = heights1[keep_ind]
    energies1_filt = energies1[keep_ind]/4.184e12
    residuals1_filt = residuals1[keep_ind]

    n = len(residuals1_filt)        

    print("Number of observations: (n = {:})".format(n))
    if R_TOL < 10:
        break

    if n < 10:
        break
    

    plt.clf()

    plt.subplot(231)

    plt.scatter((np.array(energies1[rm_ind])/4.186e12), residuals1[rm_ind], c='r')
    plt.scatter((np.array(energies1_filt)), residuals1_filt, c='k')
    plt.semilogx()
    plt.axhline(R_TOL)
    plt.xlabel("Energy [kT TNT]")
    plt.ylabel("Error in Solution [s]")

    plt.xlim([(min(energies1)/4.186e12), (max(energies1)/4.186e12)])
    plt.ylim([0, R_TOL_INIT])

    
    # fig = plt.figure()

    if CONTOUR:

        try:
            #### 2-D
            # ax = plt.subplot(232)
            # # plt.scatter(lons1_filt, lats1_filt, residuals1_filt, linewidths=0.5, colors='w')
            # cntr = plt.scatter(lons1_filt, lats1_filt, residuals1_filt, cmap="viridis")
            # plt.colorbar(cntr)

            #### 3-D
            ax = plt.subplot(232, projection="3d")
            x = lons1_filt
            y = lats1_filt
            z = heights1_filt
            cntr = ax.scatter(x, y, z, s=5, marker="o", c=residuals1_filt, cmap="viridis")
            plt.colorbar(cntr)
        except RuntimeError:
            pass

        # for stat in STATIONS:
        #     plt.scatter(stat[1], stat[0], marker="^", c="k")
        #     plt.text(stat[1], stat[0], stat[2])

    plt.xlabel("Longitude [E]")
    plt.ylabel("Latitude [N]")
    # plt.xlim([-149.5, -148.5])
    # plt.ylim([67.5, 68.5])

    plt.subplot(212)

    y_hist, x_hist, _ = plt.hist((energies1_filt), density=True, bins=int(n/BINS))
    x_hist = (x_hist[1:] + x_hist[:-1])/2 # for len(x)==len(y)


    plt.xlim([(min(energies1)/4.186e12), (max(energies1)/4.186e12)])
    # plt.ylim([0, 1])

    mu, std = norm.fit(energies1_filt)
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)

    shape, loc, scale = lognorm.fit(np.log10(energies1_filt))
    logmu = np.log(scale)
    logstd = shape

    rng = [(min(energies1_filt)), (max(energies1_filt))]

    if BIMODAL:
        try:
            params, cov = curve_fit(bimodal, x_hist, y_hist)
            plt.plot(x, bimodal(x,*params), color='red', label='Bimodal Distribution {:.2f} ({:.2f}) & {:.2f} ({:.2f}) kg TNT'\
                        .format(10**params[1], 10**params[2], 10**params[4], 10**params[5]))
        except RuntimeError:
            pass 

    #plt.plot(x, lognorm.pdf(x, shape, loc, scale), label='Log Norm: Mean {:.2f} ({:.2f}) kT TNT'.format(10**logmu, 10**logstd))
    plt.plot(x, p, 'k', linewidth=2, label="Normal Distribution: Mean {:.2E} ({:.2E}) ({:.2E} - {:.2E}) kT TNT".format(mu, std, rng[0], rng[1]))
    plt.semilogx()
    plt.title("N = {:}, R_TOL = {:}".format(n, R_TOL))
    plt.ylabel("Normalized Number of Observations")
    plt.xlabel("Energy [kT TNT]")
    plt.legend()

    print(findAvgSol(lats1_filt, lons1_filt, heights1_filt, energies1_filt, residuals1_filt))
    plt.subplot(233)

    SIGMA = 3

    r_tol_list.append(R_TOL)
    mu_list.append(mu)
    err_list.append(SIGMA*std)
    min_list.append(max([min(energies1_filt), mu - SIGMA*std]))
    max_list.append(min([max(energies1_filt), mu + SIGMA*std]))

    plt.scatter(r_tol_list, mu_list, c='k')
    plt.errorbar(r_tol_list, mu_list, yerr=err_list, fmt="o", c='k', capsize=5, label="{:} Std Dev".format(SIGMA))
    # plt.errorbar(r_tol_list, mu_list, yerr=[min_list, max_list], fmt="o", c='r', capsize=5)
    plt.semilogy()
    plt.xlabel("Allowed Tolerance [s]")
    plt.ylabel("Energy Estimate [kT TNT]")
    plt.gca().invert_xaxis()
    plt.legend()

    plt.draw()
    plt.pause(0.0001)

print(best_sol)  
plt.waitforbuttonpress()