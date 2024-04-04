#-------------------------------------------------------------------------------
# Name:        weather_datasets.py
# Purpose:     script to create weather object and other functions
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'weather_datasets.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath, join, isdir, isfile
from netCDF4 import Dataset, num2date
from glob import glob
from copy import copy
from time import sleep

WTHR_RSRC = ['CHESS']
METRICS = ['precip', 'tas', 'tasmax', 'tasmin']
REALISATIONS = list(['01', '04', '06', '15'])
RCPS = list(['rcp26', 'rcp45', 'rcp60', 'rcp85'])

RCPS = list([])
REALISATIONS = []

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '

sleepTime = 5

def _fetch_chess_wthr_nc_parms(nc_fname, wthr_rsrce, resol_time, scenario):
    """
    create object describing weather dataset characteristics
    """

    # standard names
    # ==============
    time_var_name = 'time'
    lat = 'lat'
    lon = 'lon'

    nc_fname = normpath(nc_fname)
    nc_dset = Dataset(nc_fname, 'r')
    time_var = nc_dset.variables[time_var_name]
    if 'calendar' in time_var.ncattrs():
        calendar_attr = time_var.calendar
    else:
        calendar_attr = 'standard'

    lat_var = nc_dset.variables[lat]
    lon_var = nc_dset.variables[lon]
    nrth_var = nc_dset.variables['y']
    east_var = nc_dset.variables['x']

    lats = lat_var[:][:]
    lons = lon_var[:][:]

    lat_ll = lats.min()
    lon_ll = lons.min()
    lat_ur = lats.max()
    lon_ur = lons.max()

    # resolutions
    # ===========
    resol_east = int(east_var[1] - east_var[0])
    resol_nrth = int(nrth_var[1] - nrth_var[0])
    if abs(resol_east) != abs(resol_nrth):
        print(WARN_STR + 'Check {} weather resource {}\teasting/northing resolutions: {} {} should be same'
                                                                        .format(wthr_rsrce, resol_east, resol_nrth))
    nrthngs = [int(val.item()) for val in nrth_var]
    eastngs = [int(val.item()) for val in east_var]

    # Get the start and end date of the time series (as datetime objects):
    # ====================================================================
    time_var_units = time_var.units
    start_day = int(time_var[0])
    try:
        start_date = num2date(start_day, units = time_var_units, calendar = calendar_attr)
    except (TypeError) as err:
        print(ERROR_STR + str(err) + ' deriving start and end year for dataset: ' + nc_fname)
        return None

    end_day = int(time_var[-1])
    end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)
    start_year = start_date.year
    end_year = end_date.year

    nc_dset.close()

    # construct weather resource
    # ==========================
    wthr_rsrc = {'year_start': start_year,  'year_end': end_year,
            'resol_nrth': resol_nrth, 'lat_frst': lat_ll, 'lat_last': lat_ur, 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_east': resol_east, 'lon_frst': lon_ll, 'lon_last': lon_ur, 'lon_ll': lon_ll, 'lon_ur': lon_ur,
            'nrthngs': nrthngs, 'eastngs': eastngs,
            'resol_time': resol_time,  'scenario': scenario}

    print('{} start and end year: {} {}\tresolution: {} metres'
            .format(wthr_rsrce, wthr_rsrc['year_start'],  wthr_rsrc['year_end'], abs(wthr_rsrc['resol_nrth'])))

    return wthr_rsrc

def read_wthr_dsets_detail(form, rqrd_rsces = WTHR_RSRC):
    """

    """
    wthr_rsrces_generic = list([])
    wthr_sets = {}

    form.wthr_rsrcs_generic = wthr_rsrces_generic
    form.wthr_sets = wthr_sets
    wthr_dir = form.sttngs['wthr_dir']

    # ====================
    gnrc_rsrce = 'CHESS'
    if gnrc_rsrce in rqrd_rsces:
        #  ================================= CHESS ====================================
        print()
        wthr_rsrce_base = 'CHESS_historic'

        chss_hist_flag = False
        valid_wthr_dset_rsrces = []
        chss_dir = join(wthr_dir, wthr_rsrce_base, 'Monthly')
        if isdir(chss_dir):
            chss_fns = sorted(glob(chss_dir + '/*.nc'))
            if len(chss_fns) > 0:
                wthr_rsrce = copy(wthr_rsrce_base)
                wthr_sets[wthr_rsrce] = _fetch_chess_wthr_nc_parms(chss_fns[0], wthr_rsrce, 'Monthly', 'historic')
                wthr_sets[wthr_rsrce]['base_dir'] = chss_dir
                wthr_sets[wthr_rsrce]['metrics'] = METRICS
                wthr_sets[wthr_rsrce]['fn_precip'] = glob(chss_dir + '/precip_*.nc')[0]
                wthr_sets[wthr_rsrce]['fn_tas']    = glob(chss_dir + '/tas_*.nc')[0]
                wthr_sets[wthr_rsrce]['fn_tasmax'] = glob(chss_dir + '/tasmax_*.nc')[0]
                wthr_sets[wthr_rsrce]['fn_tasmin'] = glob(chss_dir + '/tasmin_*.nc')[0]
                valid_wthr_dset_rsrces.append(wthr_rsrce)
                chss_hist_flag = True
            else:
                print('No CHESS historic datasets present in ' + chss_dir)
                sleep(sleepTime)
                exit(0)

        chss_hist_dir = copy(chss_dir)

        # check CHESS RCPs
        # ================
        chss_rcp_flag = True
        wthr_rsrce_base = 'CHESS_RCPs'
        chss_rcp_dir = join(wthr_dir, wthr_rsrce_base)
        for scenario in RCPS:
            for realis in REALISATIONS:
                chss_dir = join(chss_rcp_dir, 'Monthly', scenario + '_bias-corrected', realis)
                wthr_rsrce = 'chess_' + scenario + '_' + realis
                if isdir(chss_dir):
                    chss_fns = glob(chss_dir + '/*.nc')
                    if len(chss_fns) > 0:
                        wthr_sets[wthr_rsrce] = \
                                            _fetch_chess_wthr_nc_parms(chss_fns[0], wthr_rsrce, 'Monthly', scenario)
                        wthr_sets[wthr_rsrce]['base_dir'] = chss_dir
                        wthr_sets[wthr_rsrce]['ds_precip'] = glob(chss_dir + '/*_pr_*.nc')[0]
                        wthr_sets[wthr_rsrce]['ds_tas'] = glob(chss_dir + '/*_tas_*.nc')[0]
                        valid_wthr_dset_rsrces.append(wthr_rsrce)
                        chss_rcp_flag = True
                else:
                    print('No CHESS RCP datasets present in ' + chss_dir)

        # =========================
        if chss_hist_flag and chss_rcp_flag:
            wthr_rsrces_generic.append(gnrc_rsrce)
        else:
            if not chss_hist_flag :
                print('CHESS historic datasets incomplete in ' + chss_hist_dir)
            if not chss_rcp_flag:
                print('CHESS RCP datasets incomplete in ' + chss_rcp_dir)


    form.wthr_rsrcs_generic = wthr_rsrces_generic
    form.wthr_sets = wthr_sets

    print('')
    return

def report_aoi_size(form, lon_ll, lat_ll, lon_ur, lat_ur):
    """
    write ASCII climate files
    """
    func_name =  __prog__ + ' report_aoi_size'

    # this will be initially only NASA
    # ================================
    resource = form.combo10w.currentText()
    wthr_set = form.wthr_sets[resource]
    resol_lat = wthr_set['resol_lat']
    lat0 = wthr_set['lat0']
    resol_lon = wthr_set['resol_lon']
    lon0 = wthr_set['lon0']

    lat_indx_ll = int(round((lat_ll - lat0)/resol_lat))
    lon_indx_ll = int(round((lon_ll - lon0)/resol_lon))

    lat_indx_ur = int(round((lat_ur - lat0)/resol_lat))
    lon_indx_ur = int(round((lon_ur - lon0)/resol_lon))

    lat_indx_min = min(lat_indx_ll, lat_indx_ur)
    lat_indx_max = max(lat_indx_ll, lat_indx_ur)
    nlats = lat_indx_max - lat_indx_min + 1

    lon_indx_min = min(lon_indx_ll, lon_indx_ur)
    lon_indx_max = max(lon_indx_ll, lon_indx_ur)
    nlons = lon_indx_max - lon_indx_min + 1

    # get slice for each dataset metric
    # =================================
    mess = 'will retrieve weather for {} locations - nlats/nlons: {} x {} '.format(nlats*nlons, nlats, nlons)

    print(mess)

    return
