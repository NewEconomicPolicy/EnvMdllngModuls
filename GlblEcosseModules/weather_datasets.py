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

from cvrtcoord import WGS84toOSGB36, OSGB36toWGS84
from thornthwaite import thornthwaite

EXSTNG_WTHR_RSRCS = ['CRU', 'CHESS']
REALISATIONS = list(['01', '04', '06', '15'])
RCPS = list(['rcp26', 'rcp45', 'rcp60', 'rcp85'])
REALISATIONS = ['01', '04', '06', '15']

NGRANULARITY = 120
MONTH_NAMES_SHORT = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
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

def read_wthr_dsets_detail(form, rqrd_rsces = EXSTNG_WTHR_RSRCS):
    """
    ascertain the year span for historic datasets
    TODO: replace with approach adopted for Site Specific version of Global Ecosse
    """

    # weather set linkages
    # ====================
    wthr_rsrces_generic = list([])
    wthr_set_linkages = {}
    wthr_sets = {}

    form.wthr_rsrcs_generic = wthr_rsrces_generic
    form.wthr_set_linkages = wthr_set_linkages
    form.wthr_sets = wthr_sets
    form.wthr_settings_prev = {}

    if form.settings['weather_dir'] is None:
        return

    wthr_dir = form.settings['weather_dir']

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
                wthr_sets[wthr_rsrce]['ds_precip'] = glob(chss_dir + '/precip_*.nc')[0]
                wthr_sets[wthr_rsrce]['ds_tas'] = glob(chss_dir + '/tas_*.nc')[0]
                valid_wthr_dset_rsrces.append(wthr_rsrce)
                chss_hist_flag = True
            else:
                print('No CHESS historic datasets present in ' + chss_dir)

        chss_hist_dir = copy(chss_dir)

        # check CHESS RCPs
        # ================
        chss_rcp_flag = False
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
            wthr_set_linkages[gnrc_rsrce] = valid_wthr_dset_rsrces
        else:
            if not chss_hist_flag :
                print('CHESS historic datasets incomplete in ' + chss_hist_dir)
            if not chss_rcp_flag:
                print('CHESS RCP datasets incomplete in ' + chss_rcp_dir)

    # check CRU historic
    # ==================
    gnrc_rsrce = 'CRU'
    if gnrc_rsrce in rqrd_rsces:
        # ================================= CRU ====================================
        print()
        cru_flag = False
        valid_wthr_dset_rsrces = []
        cru_dir  = wthr_dir + '/CRU_Data'
        if isdir(cru_dir):
            wthr_rsrce = 'CRU_hist'
            cru_fnames = sorted(glob(cru_dir + '/cru*dat.nc'))
            if len(cru_fnames) > 0:
                wthr_sets[wthr_rsrce] = _fetch_wthr_nc_parms(cru_fnames[0], wthr_rsrce, 'Monthly', 'historic')
                wthr_sets[wthr_rsrce]['base_dir']   = cru_dir
                wthr_sets[wthr_rsrce]['ds_precip']  = cru_fnames[0]
                wthr_sets[wthr_rsrce]['ds_tas']     = cru_fnames[1]
                valid_wthr_dset_rsrces.append(wthr_rsrce)
                cru_flag = True
            else:
                print('No CRU historic datasets present in ' + cru_dir)

        # check ClimGen
        # =============
        climgen_flag = False
        for dset_scenario in list(['A1B','A2','B1','B2']):
            climgen_dir = join(wthr_dir, 'ClimGen', dset_scenario)
            wthr_rsrce = 'ClimGen_' + dset_scenario
            if isdir(climgen_dir):
                climgen_fnames = glob(climgen_dir + '/*.nc')
                if len(climgen_fnames) > 0:
                    wthr_sets[wthr_rsrce] = _fetch_wthr_nc_parms(climgen_fnames[0], wthr_rsrce, 'Monthly', dset_scenario)
                    wthr_sets[wthr_rsrce]['base_dir']   = climgen_dir
                    wthr_sets[wthr_rsrce]['ds_precip']  = climgen_fnames[0]
                    wthr_sets[wthr_rsrce]['ds_tas']     = climgen_fnames[1]
                    valid_wthr_dset_rsrces.append(wthr_rsrce)
                    climgen_flag = True
            else:
                print('ClimGen datasets not present in ' + climgen_dir)

        if cru_flag and climgen_flag:
            wthr_rsrces_generic.append(gnrc_rsrce)
            wthr_set_linkages[gnrc_rsrce] = valid_wthr_dset_rsrces
        else:
            if not cru_flag:
                print('CRU historic datasets incomplete in ' + cru_dir)
            if not climgen_flag:
                print('ClimGen future datasets incomplete in ' + climgen_dir)

    form.wthr_rsrcs_generic = wthr_rsrces_generic
    form.wthr_set_linkages = wthr_set_linkages
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

def write_csv_wthr_file(lgr, country, gcm_name, scenario, latitude, longitude,
                                                        start_year, end_year, pettmp_pr, pettmp_tas, out_dir):
    """
    write to file, simulation weather for the given time period
    """
    func_name =  __prog__ + ' _write_csv_wthr_file'

    metric_list = list(['Precipitation', 'Temperature','Potentional Evapotranspiration'])

    # file comprising rain and temperature
    # ====================================
    short_fname =  country + '_' + gcm_name + '_' + scenario + '.txt'
    metrics_fname = join(out_dir, short_fname)
    try:
        fhand_out = open(metrics_fname, 'w')
    except PermissionError as e:
        print(str(e))
        return

    # stanza for Potentional Evapotranspiration [mm/month]
    # ===================================================
    indx1 = 0
    pettmp_pet = []
    for year in range(start_year, end_year + 1):
        indx2 = indx1 + 12

        # temperature only is required
        # ============================
        tmean = pettmp_tas[indx1:indx2]

        # pet
        if max(tmean) > 0.0:
            pet = thornthwaite(tmean, latitude, year)
        else:
            pet = [0.0]*12
            mess = '*** Warning *** monthly temperatures are all below zero for latitude: {}\tclimate directory: {}'\
                                                                                            .format(latitude, out_dir)
            print(mess)

        pettmp_pet += [round(val,2) for val in pet]

    # identify file
    # =============
    header = 'Area of interest: {}\tLatitude: {}\tLongitude : {}\n'.format(country, latitude, longitude)
    fhand_out.write(header)

    # header for each metric
    # =======================
    header_sub = 'Month'
    for year in range(start_year, end_year + 1):
        header_sub += '\t' + str(year)

    # main loop
    # =========
    for metric, pettmp in zip(metric_list, list([pettmp_pr, pettmp_tas, pettmp_pet])):
        fhand_out.write('\n')
        fhand_out.write(metric + '\n')
        fhand_out.write(header_sub + '\n')

        recs = {}
        for mnth_nme in MONTH_NAMES_SHORT:
            recs[mnth_nme] = mnth_nme

        imnth = 0
        year = start_year
        for val in pettmp:
            mnth_nme = MONTH_NAMES_SHORT[imnth]
            recs[mnth_nme] += '\t' + str(val)
            imnth += 1
            if imnth == 12:
                imnth = 0
                year += 1
                if year > end_year:
                    break

        # write records to file
        # =====================
        for mnth_nme in MONTH_NAMES_SHORT:
            fhand_out.write(recs[mnth_nme] + '\n')

    # end of writing; close file and inform user
    # ==========================================
    nyears =  end_year - start_year + 1
    mess = 'Wrote {} years of weather data for area: {}\tto file: {}\n\tpath: {}'.format(country, nyears, short_fname,
                                                                                         out_dir)
    lgr.info(mess)
    fhand_out.close()

    return

def write_csv_wthr_file_v1(lgr, country, gcm_name, scenario, latitude, longitude,
                        start_year, end_year, pettmp_pr, pettmp_tas, out_dir):
    """
    write to file, simulation weather for the given time period
    """
    func_name =  __prog__ + ' _write_csv_wthr_file'

    # file comprising rain and temperature
    # ====================================
    short_fname =  country + '_' + gcm_name + '_' + scenario + '.txt'
    metrics_fname = join(out_dir, short_fname)
    try:
        fhand_out = open(metrics_fname, 'w')
    except PermissionError as e:
        print(str(e))
        return

    header = 'AOI\tgran_lat\tgran_lon\tlatitude\tlongitude\tyear\tmonth\tprecipitation\ttemperature\n'
    fhand_out.write(header)

    gran_lat = int(round((90.0 - latitude)*NGRANULARITY))
    gran_lon = int(round((180.0 + longitude)*NGRANULARITY))
    prefix = '{}\t{}\t{}\t{}\t{}\t'.format(country, gran_lat, gran_lon, latitude, longitude)

    # write the two metrics to file
    # =============================
    iyear = start_year
    imnth = 0
    irecs = 0
    for rain, temperature in zip(pettmp_pr, pettmp_tas):
        mnth_nme = MONTH_NAMES_SHORT[imnth]
        record = prefix + '{}\t{}\t{:.1f}\t{:.1f}\n'.format(iyear, mnth_nme, rain, temperature)
        imnth += 1
        if imnth == 12:
            imnth = 0
            iyear += 1
            if iyear > end_year:
                fhand_out.write(record)
                irecs += 1
                break

        fhand_out.write(record)
        irecs += 1

    # end of writing; close file and inform user
    # ==========================================
    mess = 'Wrote {} lines of weather data for area: {}\tto file: {}\n\tpath: {}'.format(country, irecs, short_fname,
                                                                                         out_dir)
    lgr.info(mess)
    fhand_out.close()

    return

def record_wthr_settings(scenario, hist_strt_year, hist_end_year, fut_strt_year, fut_end_year):
    """
    record weather settings
    """
    previous_settings = {'scenario': scenario, 'hist_strt_year': hist_strt_year, 'hist_end_year': hist_end_year,
                                                    'fut_strt_year': fut_strt_year, 'fut_end_year': fut_end_year}
    return previous_settings

def change_wthr_rsrc(form, wthr_rsrc = None):
    """
    during initialisation, wthr_rsrc will be specified
    otherwise get it from GUI
    """
    if wthr_rsrc == '':
        return
    if wthr_rsrc is None:
        wthr_rsrc = form.combo10w.currentText()
        '''
        form.wthr_settings_prev[wthr_rsrc] = record_wthr_settings(scenario, hist_strt_year, hist_end_year,
                                                                                        fut_strt_year, fut_end_year)
        '''

    # invoked when setting up the GUI or when there has been a change in weather resource
    # ===================================================================================
    if wthr_rsrc not in form.wthr_set_linkages:
        print('weather resource ' + wthr_rsrc + ' not in wthr_set_linkages, cannot proceed')
        return

    wthr_set_linkage = form.wthr_set_linkages[wthr_rsrc]
    wthr_set_hist = wthr_set_linkage[0]
    wthr_set_fut  = wthr_set_linkage[1]
    start_year = form.wthr_sets[wthr_set_hist]['year_start']
    end_year   = form.wthr_sets[wthr_set_hist]['year_end']
    hist_syears = list(range(start_year, end_year))
    hist_eyears = list(range(start_year + 1, end_year + 1))

    # simulation years can extend back into historical period
    # =======================================================
    start_year = min(1970, end_year)
    end_year   = form.wthr_sets[wthr_set_fut]['year_end']
    fut_syears = list(range(start_year, end_year))
    fut_eyears = list(range(start_year + 1, end_year + 1))

    # get scenarios from the future weather sets for this resource
    # ============================================================
    scenarios = []
    for wthr_set in wthr_set_linkage[1:]:
        scenarios.append(form.wthr_sets[wthr_set]['scenario'])

    # make sure scenarios are a unique list
    # ======================================
    if hasattr(form, 'combo10'):
        form.combo10.clear()
        for scenario in set(scenarios):
            form.combo10.addItem(str(scenario))
    else:
        form.combo10s.clear()
        for scenario in sorted(set(scenarios)):
            form.combo10s.addItem(str(scenario))

    form.combo10r.clear()
    for realis in REALISATIONS:
        form.combo10r.addItem(str(realis))

    form.combo09s.clear()
    for year in hist_syears:
        form.combo09s.addItem(str(year))

    form.combo09e.clear()
    for year in hist_eyears:
        form.combo09e.addItem(str(year))

    form.combo11s.clear()
    for year in fut_syears:
        form.combo11s.addItem(str(year))

    form.combo11e.clear()
    for year in fut_eyears:
        form.combo11e.addItem(str(year))

    if wthr_rsrc in form.wthr_settings_prev:
        wthr_settings_prev = form.wthr_settings_prev[wthr_rsrc]
        form.combo09s.setCurrentText(wthr_settings_prev['hist_strt_year'])
        form.combo09e.setCurrentText(wthr_settings_prev['hist_end_year'])
        if hasattr(form, 'combo10'):
            form.combo10.setCurrentText(wthr_settings_prev['scenario'])
        else:
            form.combo10s.setCurrentText(wthr_settings_prev['scenario'])
        form.combo11s.setCurrentText(wthr_settings_prev['fut_strt_year'])
        form.combo11e.setCurrentText(wthr_settings_prev['fut_end_year'])

    return

def _fetch_wthr_nc_parms(nc_fname, wthr_rsrce, resol_time, scenario):
    """
    create object describing weather dataset characteristics
    """

    # standard names
    # ==============
    time_var_name = 'time'
    if wthr_rsrce == 'NASA' or wthr_rsrce[0:5] == 'EObs_' or wthr_rsrce[0:8] == 'ClimGen_':
        lat = 'latitude'
        lon = 'longitude'
    else:
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
    lat_frst = float(lat_var[0])
    lon_frst = float(lon_var[0])
    lat_last = float(lat_var[-1])
    lon_last = float(lon_var[-1])

    if lat_last > lat_frst:
        lat_ll = lat_frst; lat_ur = lat_last
    else:
        lat_ll = lat_last; lat_ur = lat_frst

    if lon_last > lon_frst:
        lon_ll = lon_frst; lon_ur = lon_last
    else:
        lon_ll = lon_last; lon_ur = lon_frst

    # resolutions
    # ===========
    resol_lon = (lon_var[-1] - lon_var[0])/(len(lon_var) - 1)
    resol_lat = (lat_var[-1] - lat_var[0])/(len(lat_var) - 1)
    if abs(resol_lat) != abs(resol_lon):
        print('Warning - weather resource {} has different lat/lon resolutions: {} {}'
                                                        .format(wthr_rsrce, resol_lat, resol_lon))

    # Get the start and end date of the time series (as datetime objects):
    # ====================================================================
    if wthr_rsrce[0:8] == 'ClimGen_':
        # print(wthr_rsrce + ' future time units attribute: ' + time_var.units)
        start_year = int(time_var.units.split(' ')[-1])
        end_year = start_year + int(len(time_var)/12) - 1
    else:
        time_var_units = time_var.units
        start_day = int(time_var[0])
        try:
            start_date = num2date(start_day, units = time_var_units, calendar = calendar_attr)
        except (TypeError) as e:
            print('Error deriving start and end year for dataset: ' + nc_fname)
            return None

        end_day = int(time_var[-1])
        end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)
        start_year = start_date.year
        end_year = end_date.year

    longitudes = list(lon_var)
    latitudes =  []
    for lati in list(lat_var):
        latitudes.append(round(float(lati),8))  # rounding introduced for NCAR_CCSM4

    nc_dset.close()

    # construct weather resource
    # ==========================
    wthr_rsrc = {'year_start': start_year,  'year_end': end_year,
            'resol_lat': resol_lat, 'lat_frst': lat_frst, 'lat_last': lat_last, 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_lon': resol_lon, 'lon_frst': lon_frst, 'lon_last': lon_last, 'lon_ll': lon_ll, 'lon_ur': lon_ur,
            'longitudes': longitudes, 'latitudes': latitudes,
            'resol_time': resol_time,  'scenario': scenario}

    print('{} start and end year: {} {}\tresolution: {} degrees'
            .format(wthr_rsrce, wthr_rsrc['year_start'],  wthr_rsrc['year_end'], abs(wthr_rsrc['resol_lat'])))

    return wthr_rsrc