# -------------------------------------------------------------------------------
# Name:        wthr_dsets_funcs.py
# Purpose:     script to create weather object and other functions
# Author:      Mike Martin
# Created:     31/07/2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
#

__prog__ = 'wthr_dsets_funcs'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath, join, isdir, isfile
from netCDF4 import Dataset, num2date
from glob import glob

SSPS = list(['ssp126', 'ssp370'])

NGRANULARITY = 120
WARN_STR = '*** Warning *** '

def read_isimip_wthr_dsets_detail(weather_dir , gnrc_rsrce):
    """
    EFISCEN-ISIMIP
    """
    wthr_rsrces_generic = list([])
    weather_set_linkages = {}
    weather_sets = {}
    valid_wthr_dset_rsrces = []

    # check CRU historic
    # ==================
    cru_flag = False
    cru_dir = weather_dir + '/CRU_Data'
    if isdir(cru_dir):
        wthr_rsrce = 'CRU_hist'
        cru_fnames = sorted(glob(cru_dir + '/cru*dat.nc'))
        if len(cru_fnames) > 0:
            weather_sets[wthr_rsrce] = fetch_weather_nc_parms(cru_fnames[0], wthr_rsrce, 'Monthly', 'historic')
            weather_sets[wthr_rsrce]['base_dir']   = cru_dir
            weather_sets[wthr_rsrce]['ds_precip']  = cru_fnames[0]
            weather_sets[wthr_rsrce]['ds_tas']     = cru_fnames[1]
            valid_wthr_dset_rsrces.append(wthr_rsrce)
            cru_flag = True
        else:
            print('No CRU historic datasets present in ' + cru_dir)

    # check ISIMIP
    # ============
    isimip_flag = False
    for dset_scenario in SSPS:
        isimip_dir = join(weather_dir, gnrc_rsrce, 'Monthly', dset_scenario)
        wthr_rsrce = gnrc_rsrce + '_' + dset_scenario
        if isdir(isimip_dir):
            isimip_fnames = glob(isimip_dir + '/*.nc')
            if len(isimip_fnames) > 1:
                weather_sets[wthr_rsrce] = fetch_weather_nc_parms(isimip_fnames[0], wthr_rsrce, 'Monthly', dset_scenario)
                weather_sets[wthr_rsrce]['base_dir'] = isimip_dir
                weather_sets[wthr_rsrce]['ds_precip'] = isimip_fnames[1]
                weather_sets[wthr_rsrce]['ds_tas'] = isimip_fnames[0]
                valid_wthr_dset_rsrces.append(wthr_rsrce)
                isimip_flag = True
            else:
                isimip_flag = False
        else:
            isimip_flag = False

    if cru_flag and isimip_flag:
        wthr_rsrces_generic.append(gnrc_rsrce)
        weather_set_linkages[gnrc_rsrce] = valid_wthr_dset_rsrces
    else:
        if not cru_flag:
            print('CRU historic datasets incomplete in ' + cru_dir)
        if not isimip_flag:
            print('ISIMIP future datasets incomplete in ' + isimip_dir)

    print('')
    return wthr_rsrces_generic, weather_set_linkages, weather_sets

def fetch_weather_nc_parms(nc_fname, wthr_rsrce, resol_time, scenario):
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

    if wthr_rsrce.find('EObs_') == 0:
        lats = [round(float(lat), 2) for lat in list(lat_var)]  # rounding introduced for EObs
        lons = [round(float(lon), 2) for lon in list(lon_var)]
    else:
        lats = [round(float(lat), 8) for lat in list(lat_var)]  # rounding introduced for NCAR_CCSM4
        lons = [round(float(lon), 8) for lon in list(lon_var)]

    lat_frst = lats[0]
    lon_frst = lons[0]
    lat_last = lats[-1]
    lon_last = lons[-1]

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
    resol_lon = (lons[-1] - lons[0])/(len(lons) - 1)
    resol_lat = (lats[-1] - lats[0])/(len(lats) - 1)
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
        except (TypeError) as err:
            print('Error deriving start and end year for dataset: ' + nc_fname)
            return None

        end_day = int(time_var[-1])
        end_date = num2date(end_day, units = time_var_units, calendar = calendar_attr)
        start_year = start_date.year
        end_year = end_date.year

    nc_dset.close()

    # construct weather resource
    # ==========================
    wthr_rsrc = {'year_start': start_year,  'year_end': end_year,
            'resol_lat': resol_lat, 'lat_frst': lat_frst, 'lat_last': lat_last, 'lat_ll': lat_ll, 'lat_ur': lat_ur,
            'resol_lon': resol_lon, 'lon_frst': lon_frst, 'lon_last': lon_last, 'lon_ll': lon_ll, 'lon_ur': lon_ur,
            'longitudes': lons, 'latitudes': lats,
            'resol_time': resol_time,  'scenario': scenario}

    print('{} start and end year: {} {}\tresolution: {} degrees'
            .format(wthr_rsrce, wthr_rsrc['year_start'],  wthr_rsrc['year_end'], abs(wthr_rsrc['resol_lat'])))

    return wthr_rsrc