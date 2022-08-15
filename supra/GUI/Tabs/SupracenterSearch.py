import numpy as np
import matplotlib.pyplot as plt

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch
from supra.GUI.Tools.GUITools import *
from supra.Utils.AngleConv import loc2Geo
from supra.Utils.Classes import Position, Supracenter



def theoSearch(bam, prefs):
    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)
    s_info, s_name, weights = getStationData(bam.setup.station_picks_file, ref_pos)

    n_stations = len(s_name)

    xstn = s_info[0:n_stations, 0:3]

    results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, pert_num=0, theo=True)

def supSearch(bam, prefs, manual=True, results_print=False, obj=None, misfits=False, theo=False):
    """
    Function to initiate PSO Search of a Supracenter
    """
    
    if theo:
        theoSearch(bam, prefs)
        return None
    # Reference location to convert to local coordinates
    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)

    # Pull station information from picks file
    s_info, s_name, weights = getStationData(bam.setup.station_picks_file, ref_pos)

    n_stations = len(s_name)

    xstn = s_info[0:n_stations, 0:3]

    # Nominal Run
    if prefs.debug:
        print("Current status: Nominal Supracenter")

    results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=manual, pert_num=0)

    # Check for if results returns None
    try:
        # Error function is the absolute L1 norm ??
        print("Nominal Results")
        print("Error Function: {:5.2f}".format(results.f_opt))
        print("Opt: {:.4f} {:.4f} {:.2f} {:.4f}"\
            .format(results.x_opt.lat, results.x_opt.lon, results.x_opt.elev, results.motc))
    except AttributeError as e:
        errorMessage('Unable to find Supracenter Solution', 2, detail='{:}'.format(e))
        return None

    # Perturbation Runs
    if prefs.pert_en:
        pert_results = [None]*prefs.pert_num
        for i in range(prefs.pert_num):

            if prefs.debug:
                print("Current status: Perturbation {:}".format(i+1))

            pert_results[i] = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=manual, pert_num=i+1)


    
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

    if stat != 0:
        print('Residual Norm: {:.4f} s'.format(np.sqrt(norm_res)/stat))
    else:
        print('Residual Norm: Undefined')

    pert_results = []
    reses = [np.sqrt(norm_res)/stat]
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

            if stat != 0:
                print('Residual Norm: {:.4f} s'.format(np.sqrt(norm_res)/stat))
            else:
                print('Residual Norm: Undefined')
            reses.append(norm_res/stat)


    if misfits:
        genMisfits(bam, prefs, results, s_info, weights, s_name, ref_pos)


    if results_print:

        return [results, pert_results, reses, s_name]

    else:

        defTable(obj.supra_res_table, n_stations + 1, 5, headers=['Station Name', "Latitude", "Longitude", "Elevation", "Residuals"])

        if stat != 0:
            resids = norm_res/stat
        else:
            resids = np.nan

        setTableRow(obj.supra_res_table, 0, terms=["Total (Time = ({:} s))".format(results.motc), results.x_opt.lat, results.x_opt.lon, results.x_opt.elev, resids])

        for i in range(n_stations):
            
            stn_pos = Position(0, 0, 0)
            stn_pos.x, stn_pos.y, stn_pos.z = xstn[i][0], xstn[i][1], xstn[i][2]
            stn_pos.pos_geo(ref_pos)

            setTableRow(obj.supra_res_table, i + 1, terms=[s_name[i], stn_pos.lat, stn_pos.lon, stn_pos.elev, results.r[i]])

        clearLayout(obj.plots)

        supScatterPlot(bam, prefs, results, xstn, s_name, obj, manual=manual, pert_results=pert_results)

        residPlot(bam, prefs, results, pert_results, s_name, xstn, prefs.workdir, obj, manual=manual)


    return None

def normalizeResids(results, shift=0):

    norm_res = 0
    stat = 0
    for res in results.r:
        if not np.isnan(res):
            stat += 1
            norm_res += (res-shift)**2

    return norm_res/stat

def genMisfits(bam, prefs, results_nom, s_info, weights, s_name, ref_pos):
    nom_res = [results_nom.x_opt.lat, results_nom.x_opt.lon, results_nom.x_opt.elev, results_nom.motc, normalizeResids(results_nom)]

    lats = np.linspace(bam.setup.lat_min, bam.setup.lat_max)
    lons = np.linspace(bam.setup.lon_min, bam.setup.lon_max)
    elevs = np.linspace(bam.setup.elev_min, bam.setup.elev_max)
    ts = np.linspace(bam.setup.t_min, bam.setup.t_max)

    plt.figure()
    
    for ll in lats: 
        supra = Supracenter(Position(ll, results_nom.x_opt.lon, results_nom.x_opt.elev), results_nom.motc)
        results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=True, override_supra=supra)
        
        plt.scatter(ll, normalizeResids(results), color='blue')

    plt.scatter(nom_res[0], nom_res[-1], marker='*', color='red')
    plt.xlabel("Latitude [deg N]")
    plt.ylabel("Mean Station Error [s]")    

    pic_file = os.path.join(prefs.workdir, bam.setup.fireball_name, 'misfits_lat.png')
    plt.savefig(pic_file)
    plt.clf()

    for ll in lons: 
        supra = Supracenter(Position(results_nom.x_opt.lat, ll, results_nom.x_opt.elev), results_nom.motc)
        results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=True, override_supra=supra)

        plt.scatter(ll, normalizeResids(results), color='blue')
    plt.scatter(nom_res[1], nom_res[-1], marker='*', color='red')
    plt.xlabel("Longitude [deg E]")   
    plt.ylabel("Mean Station Error [s]")  
    pic_file = os.path.join(prefs.workdir, bam.setup.fireball_name, 'misfits_lon.png')
    plt.savefig(pic_file)
    plt.clf()


    for ee in elevs: 
        supra = Supracenter(Position(results_nom.x_opt.lat, results_nom.x_opt.lon, ee), results_nom.motc)
        results = psoSearch(s_info, weights, s_name, bam, prefs, ref_pos, manual=True, override_supra=supra)

        plt.scatter(ee/1000, normalizeResids(results), color='blue')
    plt.scatter(nom_res[2]/1000, nom_res[-1], marker='*', color='red')

    plt.xlabel("Elevation [km]")   
    plt.ylabel("Mean Station Error [s]")  
    pic_file = os.path.join(prefs.workdir, bam.setup.fireball_name, 'misfits_elev.png')
    plt.savefig(pic_file) 
    plt.clf()

    for tt in ts: 
        
        plt.scatter(tt, normalizeResids(results_nom, shift=tt), color='blue')
    
    plt.scatter(nom_res[3], nom_res[-1], marker='*', color='red')

    plt.xlabel("Relative Time [s]") 
    plt.ylabel("Mean Station Error [s]")  
    pic_file = os.path.join(prefs.workdir, bam.setup.fireball_name, 'misfits_time.png')
    plt.savefig(pic_file) 
    plt.clf()


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

def getStationData(picks_file, ref_pos, expectedcols=6):

    try:
        s_info, s_name, weights = convStationDat(picks_file, ref_pos)
    except TypeError as e:
        errorMessage("Unable to use station picks file!", 2, info="Error in the given file name: '{:}'".format(picks_file), detail='{:}'.format(e))
        return None, None, None
    except FileNotFoundError as e:
        errorMessage("Unable to find given station picks file!", 2, info="Could not find file: {:}".format(picks_file), detail='{:}'.format(e))
        return None, None, None
    except ValueError as e:
        errorMessage("Error in station picks file!", 2, info="Make sure that the file is formatted correctly", detail='{:}'.format(e))
        return None, None, None

    return s_info, s_name, weights

def addPlot2Parent(layout, canvas, fig):

    canvas.figure.clf()
    canvas = FigureCanvas(fig)
    canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
    layout.addWidget(canvas)    
    canvas.draw()

def supScatterPlot(bam, prefs, results, xstn, s_name, parent, manual=True, pert_results=[]):


    plt.style.use('dark_background')
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)
    ax = fig.add_subplot(1, 1, 1, projection='3d')

    ### Labels
    ax.set_title("Supracenter Locations")
    ax.set_xlabel("Latitude (deg N)", linespacing=3.1)
    ax.set_ylabel("Longitude (deg E)", linespacing=3.1)
    ax.set_zlabel('Elevation (m)', linespacing=3.1)


    n_stations = len(s_name)

    # plot station names and residuals
    for h in range(n_stations):

        # Convert station locations to geographic
        xstn[h, 0], xstn[h, 1], xstn[h, 2] = loc2Geo(bam.setup.lat_centre, bam.setup.lon_centre, 0, xstn[h, :])

        # Add station names
        ax.text(xstn[h, 0], xstn[h, 1], xstn[h, 2],  '%s' % (s_name[h]), size=10, zorder=1, color='w')


    r = results.r
    x_opt = results.x_opt
    errors = results.errors
    sup = results.sup

    c = np.nan_to_num(r)

    lat_list = []
    lon_list = []
    elev_list = []

    # Add stations with color based off of residual
    ax.scatter(xstn[:, 0], xstn[:, 1], xstn[:, 2], c=abs(c), marker='^', cmap='viridis_r', depthshade=False)

    ax.scatter(x_opt.lat, x_opt.lon, x_opt.elev, c = 'r', marker='*')
    # Try and plot trajectory if defined
    try:
        ax.plot3D([bam.setup.trajectory.pos_f.lat,  bam.setup.trajectory.pos_i.lat ],\
                  [bam.setup.trajectory.pos_f.lon,  bam.setup.trajectory.pos_i.lon ],\
                  [bam.setup.trajectory.pos_f.elev, bam.setup.trajectory.pos_i.elev],
                  'blue')
    except AttributeError:
        pass


    for ptb in range(prefs.pert_num):
        
        
        ax.scatter(pert_results[ptb].x_opt.lat, pert_results[ptb].x_opt.lon, pert_results[ptb].x_opt.elev, c = 'g', marker='*', alpha=0.7)

        lat_list.append(pert_results[ptb].x_opt.lat)   
        lon_list.append(pert_results[ptb].x_opt.lon)  
        elev_list.append(pert_results[ptb].x_opt.elev)   

    if prefs.pert_en:
        lat_list = np.array(lat_list)
        lon_list = np.array(lon_list)
        elev_list = np.array(elev_list)
        # Plot the surface
        a = np.nanmax(((x_opt.lat - lat_list)**2 + (x_opt.lon - lon_list)**2)**0.5)
        r = np.nanmax(abs(x_opt.elev - elev_list))

        x, y, z = sphereData(a, a, r, x_opt.lat, x_opt.lon, x_opt.elev)

        # uncertainty sphere
        ax.plot_wireframe(x, y, z, color='r', alpha=0.5)

    ax.text(x_opt.lat, x_opt.lon, x_opt.elev, '%s' % ('Supracenter'), zorder=1, color='w')

    if not manual:
        for i in range(len(sup)):
            sup[i, 0], sup[i, 1], sup[i, 2] = loc2Geo(bam.setup.lat_centre, bam.setup.lon_centre, 0, sup[i, :])
        sc = ax.scatter(sup[:, 0], sup[:, 1], sup[:, 2], c=errors, cmap='inferno_r', depthshade=False)
        a = plt.colorbar(sc, ax=ax)
        a.set_label("Error in Supracenter (s)")

    addPlot2Parent(parent.plots, parent.suprafig, fig)
    ax.mouse_init()


def residPlot(bam, prefs, results_arr, pert_res, s_name, xstn, output_name, parent, manual=True):
    """ outputs a 2D residual plot of the stations with the optimal supracenter

    Arguments:
        x_opt: [list] optimal supracenter position
        s_name: [list] list of station names
        xstn: [list] list of lat, lon, height positions of each station
        resid: [list] list of residuals to each station
        output_name: [string] folder to store the data in
        n_stations: [int] number of stations
    """
    n_stations = len(s_name)

    x_opt = results_arr.x_opt
    resid = results_arr.r

    # Hotfix to deal with nan stations not appearing
    new_resid = []
    for r in resid:
        if np.isnan(r):
            r = 999
        new_resid.append(r)

    resid = np.array(new_resid)

    plt.style.use('dark_background')
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)
    ax = fig.add_subplot(1, 1, 1)
    res = ax.scatter(xstn[:, 1], xstn[:, 0], c=abs(resid), marker='^', cmap='viridis_r', s=21)

    for nn in range(n_stations):
        ax.text(xstn[nn, 1], xstn[nn, 0],  '%s' % ('{:}'.format(s_name[nn])), size=10, zorder=1, color='w')


    ax.text(x_opt.lon, x_opt.lat,  '%s' % ('Supracenter'), size=10, zorder=1, color='w')

    lat_list = []
    lon_list = []

    ax.scatter(x_opt.lon, x_opt.lat, c = 'r', marker='*', s=21)
    

    if prefs.pert_en:
        for ptb in pert_res:
        
            ax.scatter(ptb.x_opt.lon, ptb.x_opt.lat, c="g", marker='*', s=21, alpha=0.7)
            lat_list.append(ptb.x_opt.lat)   
            lon_list.append(ptb.x_opt.lon)  

        lat_list = np.array(lat_list)
        lon_list = np.array(lon_list)
        # Plot the surface
        a = np.nanmax(((x_opt.lat - lat_list)**2 + (x_opt.lon - lon_list)**2)**0.5)

        circle = plt.Circle((x_opt.lon, x_opt.lat), a, color='r', alpha=0.3)
        ax.add_artist(circle)

    ax.set_ylabel("Latitude (deg N)")
    ax.set_xlabel("Longitude (deg E)")
    ax.set_title("Station Residuals")

    c = plt.colorbar(res, ax=ax)
    c.set_label("Station Residuals (s)")

    addPlot2Parent(parent.plots, parent.residfig, fig)


