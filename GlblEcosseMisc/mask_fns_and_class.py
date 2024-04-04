#-------------------------------------------------------------------------------
# Name:        mngmnt_fns_and_class.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'mngmnt_fns_and_class.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
from os.path import join, isfile, normpath
from os import remove
from netCDF4 import Dataset
from numpy import arange, float64, array, abs, unique
from math import ceil, floor
from time import strftime

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
N_DECIM = 4
HILDA_DSET = 'hilda_eu28_22_32_mask.nc'
CATEGORIES = {'0': 'Void', '11': 'Urban', '22': 'Cropland', '33': 'Pasture', '40': 'Forest (Unknown/Other)',
              '41': 'Forest (Evergreen, needle leaf)', '42': 'Forest ( Evergreen, broad leaf)',
              '43': 'Forest (Deciduous, needle leaf)', '44': 'Forest (Deciduous, broad leaf)',
              '45': 'Forest (Mixed)', '55': 'Grass/shrubland', '66': 'Other land', '77': 'Water'}

THRESHOLD = 15  # percent

def _check_voids(slices):
    '''
    check first and last years for HILDA cell
    '''
    void_flag = True
    for nyr in [0, -1]:
        slice = slices[nyr, :, :]
        uniq_res = unique(slice, return_counts=True)
        land_uses = uniq_res[0]
        if len(land_uses) == 1 and land_uses[0] == 0:
            continue
        else:
            void_flag = False
            break

    return void_flag

def fetch_prevalence(lggr, lu_defn, mask_defn, lat, lon, resol_d2, nvoids, no_squares, nerrors):
    '''
    extract data for this simulation cell
    '''
    func_name = __prog__ + ' fetch_prevalence'

    lu_totals = {key: 0 for (key, val) in CATEGORIES.items()} # dictionary comprehension


    lat_ll = lat - resol_d2
    lat_ur = lat + resol_d2
    lat_indx_min = (abs(lu_defn.lats_ar - lat_ll)).argmin()
    lat_indx_max = (abs(lu_defn.lats_ar - lat_ur)).argmin()
    if lat_indx_min > lat_indx_max:
        lat_indx = lat_indx_min
        lat_indx_min = lat_indx_max
        lat_indx_max = lat_indx

    nsize_lat = lat_indx_max - lat_indx_min

    lon_ll = lon - resol_d2
    lon_ur = lon + resol_d2
    lon_indx_min = (abs(lu_defn.lons_ar - lon_ll)).argmin()
    lon_indx_max = (abs(lu_defn.lons_ar - lon_ur)).argmin()
    nsize_lon = lon_indx_max - lon_indx_min

    if nsize_lat != nsize_lon:
        no_squares += 1

    varname = lu_defn.var_names[0]

    # only interested in last 60 years
    # ================================
    try:
        slices = lu_defn.nc_dset.variables[varname][61:, lat_indx_min:lat_indx_max, lon_indx_min:lon_indx_max]
    except RuntimeWarning as err:
        print(err)
        nerrors == 1
        return

    if _check_voids(slices):
        nvoids += 1
        return

    # check each HILDA slice from set of slices for this cell
    # =======================================================
    nslices = nsize_lat * nsize_lon
    nzeros = 0
    nbad_keys = 0
    for indx_lat in range(nsize_lat):
        for indx_lon in range(nsize_lon):

            # each slice has complete time range
            # ==================================
            slice = slices[:, indx_lat, indx_lon]
            uniq_res = unique(slice, return_counts=True)
            land_uses = uniq_res[0]
            num_lu_yrs = uniq_res[1]
            if len(land_uses) == 1 and land_uses[0] == 0:
                nzeros += 1
                continue
            else:
                out_rec = {}
                for lu, nyrs in zip(land_uses.tolist(), num_lu_yrs.tolist()):
                    land_use = str(lu)
                    out_rec[land_use] = nyrs
                    try:
                        lu_totals[land_use] += nyrs
                    except KeyError as err:
                        # print(str(err))
                        nbad_keys += 1
                        continue

                # lggr.info('\t' + str(out_rec))   # write for each HILDA slice

    # ignore uncategorised land uses ie all zero slices
    # =================================================
    if nzeros == nslices:
        nvoids += 1
    else:
        nlu_valid = (nslices - nzeros)*lu_defn.nyears - lu_totals['0']
        lu_prnct = {key: round(100 * (val / nlu_valid), 2) for (key, val) in lu_totals.items()}

        # insert into mask - ascertain correct lat/lon indices
        # ====================================================
        lat_indx, lon_indx, ret_code = mask_defn.get_nc_coords(lat, lon)
        if lu_prnct['22'] >= THRESHOLD:
            crop_land = True
        else:
            crop_land = False
        mask_defn.nc_dset.variables['cropland'][lat_indx, lon_indx] = crop_land

        if lu_prnct['33'] >= THRESHOLD:
            pasture = True
        else:
            pasture = False
        mask_defn.nc_dset.variables['pasture'][lat_indx, lon_indx] = pasture

        if lu_prnct['55'] >= THRESHOLD:
            grassland = True
        else:
            grassland = False
        mask_defn.nc_dset.variables['grassland'][lat_indx, lon_indx] = grassland

        if lu_prnct['66'] >= THRESHOLD:
            other = True
        else:
            other = False
        mask_defn.nc_dset.variables['other'][lat_indx, lon_indx] = other

        forest = False
        for ctgy in range(40,45):
            if lu_prnct[str(ctgy)] >= THRESHOLD:
                forest = True
                break

        mask_defn.nc_dset.variables['forest'][lat_indx, lon_indx] = forest

        # write summary for each cell
        # ===========================
        mess = 'Lat: {:6.2f} Lon: {:6.2f}\tN zeros: {:4d}'.format(lat, lon, nzeros)
        lggr.info(mess)

    return

def make_mask_nc(hilda_dir, resol, lon_ll, lat_ll, lon_ur, lat_ur):
    '''
    base extent on weather set and build lat long arrays
    '''
    resol_d2 = resol/2
    alons = arange(floor(lon_ll) + resol_d2, ceil(lon_ur) - resol_d2, resol, dtype=float64)
    num_alons = len(alons)
    alats = arange(floor(lat_ll) + resol_d2, ceil(lat_ur) - resol_d2, resol, dtype=float64)
    num_alats = len(alats)

    mess = 'Number of longitudes: {} and latitudes: {}'.format(num_alons, num_alats)
    print(mess)

    mask_fn = join(hilda_dir, HILDA_DSET)
    if isfile(mask_fn):
        try:
            remove(mask_fn)
        except PermissionError as err:
            print(str(err))
            return None

        print('Deleted ' + mask_fn)

    nc_dset = Dataset(mask_fn, 'w', format='NETCDF4')  # call the Dataset constructor to create file

    # create global attributes
    # ========================
    nc_dset.history = 'Land use mask for Europe consisting of Cropland and Pasture derived from the HILDA land use dataset'
    date_stamp = strftime('%H:%M %d-%m-%Y')
    nc_dset.attributation = 'Created at ' + date_stamp + ' via Spatial Ecosse'
    nc_dset.weather_dataset = HILDA_DSET
    data_used = 'Data used: HILDA land use dataset'
    nc_dset.dataUsed = data_used

    nc_dset.createDimension('latitude', num_alats)
    nc_dset.createDimension('longitude', num_alons)

    # create lat/long variables as 4 byte float
    # =========================================
    lats = nc_dset.createVariable('latitude', 'f4', ('latitude',))
    lats.description = 'degrees of latitude North to South in ' + str(resol) + ' degree steps'
    lats.units = 'degrees_north'
    lats.long_name = 'Latitude'
    lats.axis = 'Y'
    lats[:] = [round(float(lat), N_DECIM) for lat in alats]

    lons = nc_dset.createVariable('longitude', 'f4', ('longitude',))
    lons.description = 'degrees of longitude West to East in ' + str(resol) + ' degree steps'
    lons.units = 'degrees_east'
    lons.long_name = 'Longitude'
    lons.axis = 'X'
    lons[:] = [round(float(lon), N_DECIM) for lon in alons]
    
    # create cropland, pasture and forest variables
    # =============================================
    cropland = nc_dset.createVariable('cropland', 'i2', ('latitude', 'longitude'), fill_value = 0)
    cropland.units = 'none'
    cropland.missing_value = 0

    # =========================================================
    pasture = nc_dset.createVariable('pasture', 'i2', ('latitude', 'longitude'), fill_value = 0)
    pasture.units = 'none'
    pasture.missing_value = 0

    # =========================================================
    grassland = nc_dset.createVariable('grassland', 'i2', ('latitude', 'longitude'), fill_value = 0)
    grassland.units = 'none'
    grassland.missing_value = 0

    # =========================================================
    other = nc_dset.createVariable('other', 'i2', ('latitude', 'longitude'), fill_value=0)
    other.units = 'none'
    other.missing_value = 0

    # =========================================================
    forest = nc_dset.createVariable('forest', 'i2', ('latitude', 'longitude'), fill_value = 0)
    forest.units = 'none'
    forest.missing_value = 0

    # close netCDF file
    # ================
    nc_dset.sync()
    nc_dset.close()

    mess = 'Created: ' + mask_fn
    print(mess)
    # form.lgr.info()

    return mask_fn

class LanduseSet(object, ):
    '''

    '''
    def __init__(self, nc_fname):
        '''

        '''
        nc_fname = normpath(nc_fname)

        nc_dset = Dataset(nc_fname, mode='r')
        lat_var = 'latitude'
        lon_var = 'longitude'
        lats = [round(float(lat), N_DECIM) for lat in nc_dset.variables[lat_var][:]]
        lats_ar = array(lats)
        lons = [round(float(lat), N_DECIM) for lat in nc_dset.variables[lon_var][:]]
        lons_ar = array(lons)

        # record var names
        # ================
        exclude_vars = list([lat_var, lon_var, 'time'])
        year_start = None
        year_end   = None
        var_names  = []
        for var in nc_dset.variables:
            if var not in exclude_vars:
                var_names.append(var)

            if var == 'time':
                time_var = nc_dset.variables[var]
                time_units = time_var.units
                if time_units.find('years') == -1:
                    print(WARNING_STR + 'time units <' + time_units + '> not recognised in dataset ' + nc_fname)
                else:
                    year_start = int(time_var[0])
                    year_end =   int(time_var[-1])

        nc_dset.close()

        lat_frst = lats[0]
        lon_frst = lons[0]
        lat_last = lats[-1]
        lon_last = lons[-1]

        # required for bounding box
        # =========================
        if lat_last > lat_frst:
            lat_ll = lat_frst
            lat_ur = lat_last
        else:
            lat_ll = lat_last
            lat_ur = lat_frst

        if lon_last > lon_frst:
            lon_ll = lon_frst
            lon_ur = lon_last
        else:
            lon_ll = lon_last
            lon_ur = lon_frst

        self.lat_frst = lats[0]
        self.lon_frst = lons[0]
        self.lat_last = lats[-1]
        self.lon_last = lons[-1]

        self.lat_var = lat_var
        self.lon_var = lon_var
        self.bbox = lon_ll, lat_ll, lon_ur, lat_ur

        self.nc_fname = nc_fname
        self.var_names = var_names
        self.nc_dset   = None

        # resolutions
        # ===========
        self.resol_lon = round((lons[-1] - lons[0])/(len(lons) - 1), N_DECIM)
        self.resol_lat = round((lats[-1] - lats[0])/(len(lats) - 1), N_DECIM)
        self.max_lat_indx = len(lats) - 1
        self.max_lon_indx = len(lons) - 1

        #
        self.lats = lats
        self.lats_ar = lats_ar
        self.lons = lons
        self.lons_ar = lons_ar
        self.year_start = year_start
        self.year_end   = year_end
        if year_start is None:
            self.nyears = 0
        else:
            self.nyears = year_end - year_start + 1

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
