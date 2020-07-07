import numpy as np

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch

from supra.Utils.Classes import Position

def supSearch(bam, prefs):
    """
    Function to initiate PSO Search of a Supracenter
    """

    # Reference location to convert to local coordinates
    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)

    # Pull station information from picks file
    s_info, s_name, weights = getStationData(bam.setup.station_picks_file, ref_pos)

    n_stations = len(s_name)

    xstn = s_info[0:n_stations, 0:3]

    results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=True)

    print("Error Function: {:5.2f}".format(results.f_opt))
    print("Opt: {:.4f} {:.4f} {:.2f} {:.4f}"\
        .format(results.x_opt.lat, results.x_opt.lon, results.x_opt.elev, results.motc))
    
    print('Residuals:')
    for i in range(len(s_name)):
        print('{:}: {:.4f} s'.format(s_name[i], results.r[i]))

    print('Total Residual: {:.4f} s'.format(np.sum(results.r)))
    return None

    # self.scatterPlot(self.bam.setup, results, n_stations, xstn, s_name, dataset, manual=False)

    # self.residPlot(results, s_name, xstn, self.prefs.workdir, n_stations, manual=False)

    # print("Error Function: {:5.2f} (Nominal)           | Opt: {:+.4f} {:+.4f} {:.2f} {:+.4f}"\
    #     .format(results.f_opt, results.x_opt.lat, results.x_opt.lon, \
    #         results.x_opt.elev, results.motc))

    # file_name = os.path.join(self.prefs.workdir, self.bam.setup.fireball_name, "SupracenterResults.txt")
    # print('Output printed at: {:}'.format(file_name))

    # with open(file_name, "w") as f:
    #     f.write("Results\n")
    #     for ii, result in enumerate(results):
    #         if ii >= 1:
    #             f.write("Error Function: {:5.2f} (Perturbation {:4d}) | Opt: {:+.4f} {:+.4f} {:.2f} {:+.4f}\n"\
    #                 .format(results[ii].f_opt, ii, results[ii].x_opt.lat, results[ii].x_opt.lon, \
    #                     results[ii].x_opt.elev, results[ii].motc))
    #         else:
    #             f.write("Error Function: {:5.2f} (Nominal)           | Opt: {:+.4f} {:+.4f} {:.2f} {:+.4f}\n"\
    #                 .format(results[ii].f_opt, results[ii].x_opt.lat, results[ii].x_opt.lon, \
    #                     results[ii].x_opt.elev, results[ii].motc))

    # defTable(self.sup_results_table, n_stations + 1, 5, headers=['Station Name', "Latitude", "Longitude", "Elevation", "Residuals"])

    # setTableRow(self.sup_results_table, 0, terms=["Total (Time = ({:}s)".format(results[0].motc), results[0].x_opt.lat, results[0].x_opt.lon, results[0].x_opt.elev, results[0].f_opt])

    # for i in range(n_stations):
    #     setTableRow(self.sup_results_table, i + 1, terms=[s_name[i], xstn[i][0], xstn[i][1], xstn[i][2], results[0].r[i]])

def getStationData(picks_file, ref_pos):

    try:
        s_info, s_name, weights = convStationDat(picks_file, ref_pos)
    except TypeError as e:
        errorMessage("Unable to use station picks file!", 2, info="Error in the given file name: '{:}'".format(self.bam.setup.station_picks_file), detail='{:}'.format(e))
        return None, None, None
    except FileNotFoundError as e:
        errorMessage("Unable to find given station picks file!", 2, info="Could not find file: {:}".format(self.bam.setup.station_picks_file), detail='{:}'.format(e))
        return None, None, None
    except ValueError as e:
        errorMessage("Error in station picks file!", 2, info="Make sure that the file is formatted correctly", detail='{:}'.format(e))
        return None, None, None

    return s_info, s_name, weights