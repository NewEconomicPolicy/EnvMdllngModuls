#-------------------------------------------------------------------------------
# Name:        country_fns.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'country_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
# 
import os
import sys
from time import time
from locale import format_string
from netCDF4 import Dataset
from numpy import arange, float64, array, abs, unique, ndarray, where
from math import ceil, floor
from pandas import read_csv
from time import strftime

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
N_DECIM = 4

def update_progress(last_time, start_time, completed, est_num_sims, skipped, warning_count):

    """Update progress bar."""
    new_time = time()
    if new_time - last_time > 5:
        # display commas
        # ==============
        remain_str = format_string("%10d", est_num_sims - completed, grouping=True)
        mess = '\rCompleted: {:=6d} Skipped: {:=5d} Warnings: {:=5d} Remaining: {}'\
                                    .format(completed, skipped, warning_count, remain_str)
        sys.stdout.flush()
        sys.stdout.write(mess)
        last_time = new_time

    return last_time

def make_country_nc(csv_fname, states_fn, resol):
    '''
    base extent on weather set and build lat long arrays
    '''
    print('Creating dataframe from ' + csv_fname)
    dummy, orig_dset = os.path.split(csv_fname)
    headers = ['lon', 'lat', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9', 'cntry_code']
    data_frame = read_csv(csv_fname, skiprows=1, names=headers, delim_whitespace=True)

    alats = unique(data_frame['lat'])
    num_alats = len(alats)
    alons = unique(data_frame['lon'])
    num_alons = len(alons)

    mess = 'Number of longitudes: {} and latitudes: {}'.format(num_alons, num_alats)
    print(mess)

    nc_dset = Dataset(states_fn, 'w', format='NETCDF4')  # call the Dataset constructor to create file

    # create global attributes
    # ========================
    nc_dset.history = 'World countries'
    date_stamp = strftime('%H:%M %d-%m-%Y')
    nc_dset.attributation = 'Created at ' + date_stamp + ' via Spatial Ecosse'
    nc_dset.original_dataset = orig_dset
    nc_dset.resolution = resol
    nc_dset.dataUsed = 'Data used: from Astley Hastings supplied by Shell'

    nc_dset.createDimension('latitude', num_alats)
    nc_dset.createDimension('longitude', num_alons)

    # create lat/long variables as 4 byte float
    # =========================================
    lats = nc_dset.createVariable('latitude', 'f4', ('latitude',))
    lats.description = 'degrees of latitude South to North in ' + str(resol) + ' degree steps'
    lats.units = 'degrees_north'
    lats.long_name = 'Latitude'
    lats.axis = 'Y'
    lats[:] = alats

    lons = nc_dset.createVariable('longitude', 'f4', ('longitude',))
    lons.description = 'degrees of longitude West to East in ' + str(resol) + ' degree steps'
    lons.units = 'degrees_east'
    lons.long_name = 'Longitude'
    lons.axis = 'X'
    lons[:] = alons
    
    # create countries variable
    # =========================
    countries = nc_dset.createVariable('countries', 'i2', ('latitude', 'longitude'), fill_value = 0)
    countries.units = 'none'
    countries.missing_value = 0

    # cycle through all records
    # =========================
    print('Populating countries variable from dataframe...')
    nrecs = len(data_frame['lat'])
    skipped = 0
    warning_count = 0
    last_time = time()
    start_time = time()
    completed = 0
    for lat, lon, cntry_code in zip(data_frame['lat'], data_frame['lon'], data_frame['cntry_code']):
        res = where(alats == lat)
        lat_indx = res[0][0]
        res = where(alons == lon)
        lon_indx = res[0][0]
        countries[lat_indx, lon_indx] = int(cntry_code)
        completed += 1
        last_time = update_progress(last_time, start_time, completed, nrecs, skipped, warning_count)

    # close netCDF file
    # ================
    nc_dset.sync()
    nc_dset.close()

    mess = 'Created: ' + states_fn
    print(mess)
    # form.lgr.info()

    return
