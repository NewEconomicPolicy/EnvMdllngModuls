#-------------------------------------------------------------------------------
# Name:        filter_hwsd_fns
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'filter_hwsd_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from glob import glob
from os.path import join, isfile

from mngmnt_fns_and_class import ManagementSet
from netCDF4 import Dataset
from country_fns import update_progress
from time import time

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
N_DECIM = 4

def expand_phenology(form):
    '''
    read CSV file and add this phenology to existing NC file
    '''
    crop_mapped = 'Wheat_Spring'

    # NC dataset
    # ==========
    pheno_path = join(form.proj_path, 'Phenology')
    fnames = glob(pheno_path + '\\' + crop_mapped + '*.nc')
    if len(fnames) == 0:
        print(ERROR_STR + 'NC file for ' + crop_mapped + ' must exist')
        return

    pheno_nc_fname = fnames[0]

    # CSV dataset
    # ===========
    phnlgy_csv = join(pheno_path, 'vault','SOWING_HARVEST_DATES_SWHEAT_FILTERED.txt')
    if not isfile(phnlgy_csv):
        print(ERROR_STR + 'file ' + phnlgy_csv + ' must exist')
        return

    print('Adding phenology records from ' + phnlgy_csv + '\n\tto ' + pheno_nc_fname)

    resource = 'phenology'
    pheno_defn = ManagementSet(pheno_nc_fname, resource)
    pheno_defn.nc_dset = Dataset(pheno_defn.nc_fname, mode='a')
    pheno_defn.nc_dset.close()

    # =========
    sowing_crop_type = 'spring wheat'
    form.crop_type_sowing = {'type': sowing_crop_type, 'frst_year': -999, 'last_year': -999}
    # data_frame = read_sowing_dates_file_short(form, phnlgy_csv)

    return

def filter_hwsd_csv(form):
    '''
    base extent on weather set and build lat long arrays
    '''
    mask_fn = 'E:\\GlobalEcosseData\\Hilda_land_use\\hilda_eu28_22_32_mask.nc'
    resource = 'cropmasks'
    varname = 'cropland'
    mask_defn = ManagementSet(mask_fn, resource)
    mask_defn.nc_dset = Dataset(mask_defn.nc_fname, mode='r')

    nzeroes = 0
    n_ones = 0

    last_time = time()
    start_time = time()
    completed = 0
    df_len = len(form.hwsd_mu_globals.data_frame)

    for rec in form.hwsd_mu_globals.data_frame.values:

        lat = rec[3]
        lon = rec[4]
        lat_indx, lon_indx, ret_code = mask_defn.get_nc_coords(lat, lon)
        val = mask_defn.nc_dset.variables[varname][lat_indx, lon_indx]
        ival = int(val.item())
        if ival == 0:
            nzeroes += 1
        else:
            n_ones += 1

        completed += 1

        last_time = update_progress(last_time, start_time, completed, df_len, nzeroes, n_ones)

    mask_defn.nc_dset.close()

    mess = '\nFound - zeroes: {}\tones: {}\ttotal: {}\tdf total: {}'.format(nzeroes, n_ones, nzeroes + n_ones, df_len)
    print(mess)

    return
