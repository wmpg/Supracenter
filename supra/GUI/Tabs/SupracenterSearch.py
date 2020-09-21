import numpy as np

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch

from supra.Utils.Classes import Position

def supSearch(bam, prefs, manual=True, results=False):
    """
    Function to initiate PSO Search of a Supracenter
    """
    
    # Reference location to convert to local coordinates
    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)

    # Pull station information from picks file
    s_info, s_name, weights = getStationData(bam.setup.station_picks_file, ref_pos)

    n_stations = len(s_name)

    xstn = s_info[0:n_stations, 0:3]

    # TODO Test this code
    # TODO Display results in a statistically significant way

    # Nominal Run
    if prefs.debug:
        print("Current status: Nominal Supracenter")

    results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=manual, pert_num=0)

    # Perturbation Runs
    if prefs.pert_en:
        pert_results = [None]*prefs.pert_num
        for i in range(prefs.pert_num):

            if prefs.debug:
                print("Current status: Perturbation {:}".format(i+1))

            pert_results[i] = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=manual, pert_num=i+1)

    # Error function is the absolute L1 norm ??
    print("Nominal Results")
    print("Error Function: {:5.2f}".format(results.f_opt))
    print("Opt: {:.4f} {:.4f} {:.2f} {:.4f}"\
        .format(results.x_opt.lat, results.x_opt.lon, results.x_opt.elev, results.motc))
    
    # Results show an L2 norm, normalized by the number of stations that have arrivals
    print('Residuals:')
    for i in range(len(s_name)):
        print('{:}: {:.4f} s'.format(s_name[i], results.r[i]))

    norm_res = 0
    stat = 0
    for res in results.r:
        if not np.isnan(res):
            stat += 1
            norm_res += res**2

    print('Residual Norm: {:.4f} s'.format(norm_res/stat))
    reses = [norm_res/stat]
    if prefs.pert_en:
        for i in range(prefs.pert_num):

            # Error function is the absolute L1 norm ??
            print("Perturbation {:} Results".format(i+1))
            print("Error Function: {:5.2f}".format(pert_results[i].f_opt))
            print("Opt: {:.4f} {:.4f} {:.2f} {:.4f}"\
                .format(pert_results[i].x_opt.lat, pert_results[i].x_opt.lon, pert_results[i].x_opt.elev, pert_results[i].motc))
            
            # Results show an L2 norm, normalized by the number of stations that have arrivals
            print('Residuals:')
            for ii in range(len(s_name)):
                print('{:}: {:.4f} s'.format(s_name[ii], pert_results[i].r[ii]))

            norm_res = 0
            stat = 0
            for res in pert_results[i].r:
                if not np.isnan(res):
                    stat += 1
                    norm_res += res**2

            print('Residual Norm: {:.4f} s'.format(norm_res/stat))
            reses.append(norm_res/stat)


    if results:
        return [results, pert_results, reses, s_name]

    return None

def resultsPrint(results, pert_results, reses, s_name, prefs, doc=None):
    
    print("#####################################")
    print("SUMMARY")
    print("#####################################")

    lats = []
    lons = []
    elevs = []
    times = []

    for i in range(prefs.pert_num):
        lats.append(pert_results[i].x_opt.lat)
        lons.append(pert_results[i].x_opt.lon)
        elevs.append(pert_results[i].x_opt.elev)
        times.append(pert_results[i].motc)

    max_lat = np.nanmax(lats) - results.x_opt.lat
    max_lon = np.nanmax(lons) - results.x_opt.lon
    max_elev = np.nanmax(elevs) - results.x_opt.elev
    max_time = np.nanmax(times) - results.motc

    min_lat = results.x_opt.lat - np.nanmin(lats)
    min_lon = results.x_opt.lon - np.nanmin(lons)
    min_elev = results.x_opt.elev - np.nanmin(elevs)
    min_time = results.motc - np.nanmin(times)
    print('')
    print("Latitude : {:.4f} +{:.4f} / -{:.4f} deg N".format(results.x_opt.lat, max_lat, min_lat))
    print("Longitude: {:.4f} +{:.4f} / -{:.4f} deg E".format(results.x_opt.lon, max_lon, min_lon))
    print("Elevation: {:.4f} +{:.4f} / -{:.4f} km".format(results.x_opt.elev/1000, max_elev/1000, min_elev/1000))
    print("Time     : {:.4f} +{:.4f} / -{:.4f} s".format(results.motc, max_time, min_time))
    print('')
    print('Residuals')
    print('------------------------------------')
    print("Total Normed Residuals: {:.4f} +{:.4f} / -{:.4f} s".format(reses[0], np.nanmax(reses), np.nanmin(reses)))
    print("By Station (Nominal):")
    for ii in range(len(s_name)):
        print('{:}: {:.4f} s'.format(s_name[ii], results.r[ii]))

    if doc is not None:
        doc.add_heading('Supracenter', 2)

        doc.add_paragraph("Latitude : {:.4f} +{:.4f} / -{:.4f} deg N".format(results.x_opt.lat, max_lat, min_lat))
        doc.add_paragraph("Longitude: {:.4f} +{:.4f} / -{:.4f} deg E".format(results.x_opt.lon, max_lon, min_lon))
        doc.add_paragraph("Elevation: {:.4f} +{:.4f} / -{:.4f} km".format(results.x_opt.elev/1000, max_elev/1000, min_elev/1000))
        doc.add_paragraph("Time     : {:.4f} +{:.4f} / -{:.4f} s".format(results.motc, max_time, min_time))
        doc.add_paragraph('')
        doc.add_paragraph('Residuals')
        doc.add_paragraph("Total Normed Residuals: {:.4f} +{:.4f} / -{:.4f} s".format(reses[0], np.nanmax(reses), np.nanmin(reses)))
        doc.add_paragraph("By Station (Nominal):")
        for ii in range(len(s_name)):
            doc.add_paragraph('{:}: {:.4f} s'.format(s_name[ii], results.r[ii]))

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