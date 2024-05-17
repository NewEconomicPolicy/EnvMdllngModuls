#-------------------------------------------------------------------------------
# Name:        spec_NCfuncs.py
# Purpose:     Functions to create and write to netCDF files and return latitude and longitude indices
# Author:      Mike Martin
# Created:     25/01/2017
# Description: create dimensions: "longitude", "latitude" and "time"
#              create five ECOSSE variables i.e. 'n2o','soc','co2', 'no3', and 'ch4'
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'netcdf_funcs.py'
__version__ = '0.0.0'
__author__ = 's03mm5'

from os.path import normpath, isfile, join, split
from os import remove
import time
from netCDF4 import Dataset
from numpy import arange, float32

missing_value = -999.0
imiss_value = int(missing_value)
WARNING_STR = '*** Warning *** '

def get_nc_coords(bbox_nc, latitude, longitude, max_lat_indx, max_lon_indx):
    """

    """
    resol = 0.5
    ll_lon, ll_lat, ur_lon, ur_lat = bbox_nc
    lat_indx = round((latitude  - ll_lat)/resol)
    lon_indx = round((longitude - ll_lon)/resol)

    if lat_indx < 0 or lat_indx > max_lat_indx:
        print()
        print(WARNING_STR + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'
                                                            .format(lat_indx, round(latitude, 4), max_lat_indx))
        return -1, -1

    if lon_indx < 0 or lon_indx > max_lon_indx:
        print()
        print(WARNING_STR + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'
                                                            .format(lon_indx, round(longitude, 4), max_lon_indx))
        return -1, -1

    return lat_indx, lon_indx

def create_soil_nc_dset(form, ll_lon, ll_lat, ur_lon, ur_lat, soil_metrics):
    """

    """
    func_name =  __prog__ + ' create_soil_nc_dset'

    delete_flag = True

    # expand bounding box to make sure all results are included
    # =========================================================   
    resol = 0.5
    bbox = list([ll_lon, ll_lat, ur_lon, ur_lat])

    # build lat long arrays
    # =====================
    alons = arange(ll_lon, ur_lon + resol, resol, dtype=float32)
    alats = arange(ll_lat, ur_lat + resol, resol, dtype=float32)
    num_alons = len(alons)
    num_alats = len(alats)
    nlayers = 2
    form.bbox_nc = bbox

    mess = 'Number of longitudes: {} and latitudes: {}\n'.format(num_alons, num_alats)
    print(mess)       

    # construct the output file name and delete if it already exists
    # ==============================================================
    soil_dir = join(split(form.sims_dir)[0], 'soil_metrics')
    fout_name = normpath(join(soil_dir, 'soil_vars.nc'))
    if isfile(fout_name):
        if delete_flag:
            try:
                remove(fout_name)
                print('Deleted file: ' + fout_name)
            except PermissionError:
                print('Function: {}\tcould not delete file: {}'.format(func_name, fout_name))
                return 1
        else:
            return fout_name

    # call the Dataset constructor to create file
    # ===========================================
    nc_dset = Dataset(fout_name,'w', format='NETCDF4')

    # create global attributes
    # ========================
    history = 'Soil details for AOI extent - longitude: {} to {}'.format(ll_lon, ur_lon)
    history += '\tlatitude: {} to {}'.format(ll_lat, ur_lat)
    nc_dset.history = history
    date_stamp = time.strftime('%H:%M %d-%m-%Y')
    nc_dset.attributation = 'Created at ' + date_stamp + ' from Spatial Ecosse '
    data_used = 'Data used: HWSD soil'
    nc_dset.dataUsed = data_used

    # we create the dimensions using the createDimension method of a Group (or Dataset) instance
    nc_dset.createDimension('lat', num_alats)
    nc_dset.createDimension('lon', num_alons)
    nc_dset.createDimension('layer', nlayers)

    '''
    create the variable (4 byte float in this case)
    to create a netCDF variable, use the createVariable method of a Dataset (or Group) instance.
    first argument is name of the variable, second is datatype, third is a tuple with the name (s) of the dimension(s).
    e.g. lats = nc_dset.createVariable('latitude',dtype('float32').char,('lat',))
    '''
    lats = nc_dset.createVariable('latitude','f4',('lat',))
    lats.units = 'degrees of latitude North to South in ' + str(resol) + ' degree steps'
    lats.long_name = 'latitude'
    lats[:] = alats

    lons = nc_dset.createVariable('longitude','f4',('lon',))
    lons.units = 'degrees of longitude West to East in ' + str(resol) + ' degree steps'
    lons.long_name = 'longitude'
    lons[:] = alons

    #
    layers = nc_dset.createVariable('layer','i2',('layer',))
    layers.units = 'layer 1 is topsoil 0 to 30 cms, layer 2 is subsoil 30 to 100 cms'
    layers[:] = [1,2]

    #
    mu_globals = nc_dset.createVariable('mu_global', 'i2', ('lat','lon'), fill_value = missing_value)
    mu_globals.units = 'HWSD mapping unit'

    # create the soil variables
    # =========================
    for metric in soil_metrics:
        var_varia = nc_dset.createVariable(metric,'f4',('lat','lon','layer'), fill_value = missing_value)
        var_varia.units = soil_metrics[metric]

    # close netCDF file
    # ================
    nc_dset.sync()
    nc_dset.close()

    mess = 'Created {} netCDF file'.format(fout_name)
    print(mess)

    return fout_name
