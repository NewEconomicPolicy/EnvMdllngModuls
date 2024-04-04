#-------------------------------------------------------------------------------
# Name:        mngmnt_fns_and_class.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'dset_fns_and_class.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import normpath
from netCDF4 import Dataset
from numpy import ma
from calendar import monthrange

ERROR_STR = '*** Error *** '
ELEV_VAR_NAMES = list(['lat', 'lon', 'Band1'])
GRANULARITY = 120
numSecsDay = 3600*24

def fetch_wthr_elev(fobj, easting, nrthing, elev_defn, wthr_defn):
    """
    get precipitation or temperature data for a given variable and lat/long index for all times
    """
    func_name = __prog__ + ' fetch_chess_NC_data'

    strt_yr = wthr_defn['year_start']

    pettmp = {}
    elev, lat, lon = 3*[None]

    indx_nrth = nrthing/1000
    indx_east = easting/1000

    for metric in wthr_defn['metrics']:
        nc_dset = wthr_defn['ds_' + metric]
        try:
            if strt_yr == 1980:
                slice = nc_dset.variables[metric][:-1, indx_nrth, indx_east]
            else:
                slice = nc_dset.variables[metric][:, indx_nrth, indx_east]

        except (RuntimeWarning, KeyError) as err:
            print(ERROR_STR + err + ' at E/N: {} {}\tfor metric: '.format(easting, nrthing, metric))
            pettmp = None
            break

        if ma.is_masked(slice):
            slice_is_masked_flag = True
            fobj.write('Slice is masked at E/N: {} {}\tfor metric: '.format(easting, nrthing, metric))
            pettmp = None
            break

        # generate days per month
        # ======================
        if metric == 'precip':
            lat = float(nc_dset.variables['lat'][indx_nrth, indx_east])
            lon = float(nc_dset.variables['lon'][indx_nrth, indx_east])

            lat_indx = round((lat - elev_defn.lat_frst) / elev_defn.resol_lat)
            lon_indx = round((lon - elev_defn.lon_frst) / elev_defn.resol_lon)
            elev = float(elev_defn.nc_dset.variables['Band1'][lat_indx, lon_indx])

            days_per_month = []
            nmonths = len(slice)
            for year in range(strt_yr, strt_yr + int(nmonths / 12)):
                for imnth in range(12):
                    dummy, ndays = monthrange(year, imnth + 1)
                    days_per_month.append(ndays)

        # reform slice
        # ============
        lat = float(nc_dset.variables['lat'][indx_nrth, indx_east])
        gran_lat = round((90.0 - lat) * GRANULARITY)

        lon = float(nc_dset.variables['lon'][indx_nrth, indx_east])
        gran_lon = round((180.0 + lon) * GRANULARITY)

        if metric == 'precip':

            # convert from units = kg m-2 s-1 to mm
            # =====================================
            precip_mm = []
            for ndays, precip in zip(days_per_month, slice):
                val_mm = float(precip) * numSecsDay * ndays
                precip_mm.append(val_mm)

            pettmp[metric] = [round(val, 3) for val in precip_mm]
        else:
            pettmp[metric] = [round(val - 273.15, 3) for val in slice]

    return pettmp, elev, lat, lon

def open_proj_NC_sets(elev_defn, wthr_defn):
    """
    """
    elev_defn.nc_dset = Dataset(elev_defn.nc_fname, mode='r')

    for metric in wthr_defn['metrics']:
        wthr_defn['ds_' + metric] = Dataset(wthr_defn['fn_' + metric], mode='r')

    return

def close_proj_NC_sets(elev_defn, wthr_defn):
    """
    """
    elev_defn.nc_dset.close()
    for metric in wthr_defn['metrics']:
        wthr_defn['ds_' + metric].close()

class ElevationSet(object, ):
    """

    """
    def __init__(self, nc_fname):
        """

        """
        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname, mode='r')

        for var_name in ELEV_VAR_NAMES:

            if var_name not in nc_dset.variables:
                print(ERROR_STR + 'Bad dataset ' + nc_fname)
                self.nc_dset   = None
        else:
            lat_var = 'lat'
            lon_var = 'lon'
            lats = nc_dset.variables[lat_var][:]
            lons = nc_dset.variables[lon_var][:]

            var_names  = ['Band1']
            nc_dset.close()

            nlons = len(lons)
            nlats = len(lats)

            lat_ll = float(lats[0])
            lon_ll = float(lons[0])
            lat_ur = float(lats[-1])
            lon_ur = float(lons[-1])

            self.lat_var = lat_var
            self.lon_var = lon_var
            self.bbox =  lon_ll, lat_ll, lon_ur, lat_ur

            self.nc_fname = nc_fname
            self.var_names = var_names

            # resolutions
            # ===========
            self.resol_lon = (lons[-1] - lons[0])/(nlons - 1)
            self.resol_lat = (lats[-1] - lats[0])/(nlats - 1)
            self.max_lat_indx = nlats - 1
            self.max_lon_indx = nlons - 1

            #
            self.lats = list(lats)
            self.lons = list(lons)

            self.lat_frst = float(lats[0])
            self.lon_frst = float(lons[0])
            self.lat_last = float(lats[-1])
            self.lon_last = float(lons[-1])

def get_nc_coords(self, latitude, longitude):

        ret_code = 'OK'

        lat_indx = int(round((latitude -  self.lat_frst)/self.resol_lat))
        lon_indx = int(round((longitude - self.lon_frst)/self.resol_lon))

        if lat_indx < 0 or lat_indx > self.max_lat_indx:
            ret_code = '*** Warning *** latitude index {} out of bounds for latitude {}\tmax indx: {}'.format(lat_indx,
                                                                                round(latitude, 4), self.max_lat_indx)
        if lon_indx < 0 or lon_indx > self.max_lon_indx:
            ret_code = '*** Warning *** longitude index {} out of bounds for longitude {}\tmax indx: {}'.format(lon_indx,
                                                                                round(longitude, 4), self.max_lon_indx)
        return lat_indx, lon_indx, ret_code
