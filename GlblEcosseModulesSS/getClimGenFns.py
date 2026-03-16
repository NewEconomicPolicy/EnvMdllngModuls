"""
#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
# -------------------------------------------------------------------------------
"""
__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

from numpy.ma.core import MaskedConstant, MaskError
from netCDF4 import Dataset
from warnings import filterwarnings
from time import time
from sys import stdout

ERROR_STR = '*** Error *** '
WARNING = '*** Warning *** '

NULL_VALUE = -9999
GRANULARITY = 120

def compress_pettmp_hist_fut(pettmp_hist, pettmp_fut):
    """
    compress historic and future weather
    """
    func_name = 'compress_pettmp_hist_fut'
    keys_hist = list(pettmp_hist['precip'].keys())
    keys_fut = list(pettmp_fut['precip'].keys())

    mess = '\nHistoric and Future pettmp have '
    if keys_fut.sort() == keys_hist.sort():
        all_keys = keys_fut
        print(mess + 'same coordinate set')
    else:
        in_fut = set(keys_fut)
        in_hist = set(keys_hist)
        in_fut_not_in_hist = in_fut - in_hist
        all_keys = in_hist + list(in_fut_not_in_hist)
        print(WARNING + mess + 'same different coordinate set')

    # main loop
    # =========
    pettmp_hist_new = {'lat_lons': {}, 'precip': {}, 'tas': {}}
    pettmp_fut_new = {'precip': {}, 'tas': {}}
    nskipped = 0
    ncopied = 0
    for key in all_keys:
        if key in keys_fut and key in keys_hist:
            if pettmp_hist['precip'][key] is None:
                nskipped += 1
                continue

            pettmp_hist_new['lat_lons'][key] = pettmp_hist['lat_lons'][key]
            for metric in ['precip', 'tas']:
                pettmp_hist_new[metric][key] = pettmp_hist[metric][key]
                pettmp_fut_new[metric][key] = pettmp_fut[metric][key]
            ncopied += 1

    print('Copied {} and skipped {} weather cells in function {}'.format(ncopied, nskipped, func_name))

    return pettmp_hist_new, pettmp_fut_new

def fetch_WrldClim_NC_data(lgr, aoi_indices, climgen, nc_dsets, hist_flag=False, fut_start_indx=0, report_flag=False):
    """
    self, aoi_indices, num_band, fut_start_indx=0
    get precipitation or temperature data for a given variable and lat/long indices for all times
    """
    func_name = __prog__ + ' fetch_WrldClim_NC_data'
    filterwarnings("error")

    nkey_masked = 0
    lat_indx_min, lat_indx_max, lon_indx_min, lon_indx_max = aoi_indices
    ncells = (lat_indx_max + 1 - lat_indx_min) * (lon_indx_max + 1 - lon_indx_min)
    pettmp = {}
    pettmp['lat_lons'] = {}
    last_time = time()

    for metric in list(['precip', 'tas']):
        pettmp[metric] = {}

        if hist_flag:
            varname = climgen.hist_wthr_set_defn[metric]
            try:
                slice = nc_dsets[metric].variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]
            except BaseException as err:
                print(ERROR_STR + str(err))
                return None
        else:
            varname = climgen.fut_wthr_set_defn[metric]
            slice = nc_dsets[metric].variables[varname][:, lat_indx_min:lat_indx_max + 1, lon_indx_min:lon_indx_max + 1]

        # reform slice
        # ============
        icells = 0
        for ilat, lat_indx in enumerate(range(lat_indx_min, lat_indx_max + 1)):
            lat =  float(nc_dsets[metric].variables['lat'][lat_indx])
            gran_lat = round((90.0 - lat) * GRANULARITY)

            for ilon, lon_indx in enumerate(range(lon_indx_min, lon_indx_max + 1)):
                try:
                    lon = float(nc_dsets[metric].variables['lon'][lon_indx])
                except IndexError as err:
                    continue

                gran_lon = round((180.0 + lon) * GRANULARITY)
                key = '{:0=5d}_{:0=5d}'.format(int(gran_lat), int(gran_lon))
                icells += 1
                if report_flag:
                    last_time = update_fetch_progress(last_time, nkey_masked, icells, ncells)

                # validate values
                # ===============
                pettmp[metric][key] = NULL_VALUE
                val = slice[0, ilat, ilon]
                if type(val) is MaskedConstant:
                    lgr.info('val is ma.masked for key ' + key)
                    pettmp[metric][key] = None
                    nkey_masked += 1

                # add data for this coordinate
                # ============================
                if pettmp[metric][key] == NULL_VALUE:
                    record = [round(val, 1) for val in slice[:, ilat, ilon]]
                    pettmp[metric][key] = record[fut_start_indx:]

                pettmp['lat_lons'][key] = [lat, lon]

    return pettmp

def _apply_start_year_correction(sim_strt_yr, hist_dset_defn, pettmp):
    """
    assume monthly datasets
    check and, if necessary, correct situation where sim_strt_yr is before historic dataset start year
    """
    repeat_period = (hist_dset_defn['year_start'] - sim_strt_yr)*12
    if repeat_period <= 0:
        return pettmp

    new_pettmp = {}
    for metric in pettmp:
        new_pettmp[metric] = pettmp[metric][0:repeat_period] + pettmp[metric]

    return new_pettmp

def _fetch_wthrset_indices(wthr_set_defn, sim_strt_yr, sim_end_yr):
    """
    get indices for simulation years for monthly weather set
    """
    wthr_yr_strt = wthr_set_defn['year_start']
    wthr_yr_end = wthr_set_defn['year_end']

    # simulation end year is before start of this dataset - nothing to do
    # ===================================================================
    if sim_end_yr < wthr_yr_strt:
        return 3*[None]

    # simulation start year is after end of this dataset - nothing to do
    # ===================================================================
    if wthr_yr_strt > sim_end_yr:
        return 3 * [None]

    indx_strt = max(0, (sim_strt_yr - wthr_yr_strt)*12)

    if sim_end_yr >= wthr_yr_end:

        # simulation end year is in future and beyond this dataset end year
        # =================================================================
        indx_end = -1
        next_strt_yr = wthr_yr_end + 1
    else:
        # simulation end year is before this dataset end year
        # ===================================================
        indx_end = (sim_end_yr - wthr_yr_end)*12
        next_strt_yr = -1

    return indx_strt, indx_end, next_strt_yr

def join_hist_fut_to_sim_wthr(climgen, pettmp_hist, pettmp_fut):
    """
    join historic and future weather
    TODO: can be made more efficient by doing this once
    """
    sim_strt_yr = climgen.sim_start_year
    sim_end_yr = climgen.sim_end_year
    indx_hist_strt, indx_hist_end, next_strt_yr = _fetch_wthrset_indices(climgen.hist_wthr_set_defn,
                                                                            sim_strt_yr, sim_end_yr)
    indx_fut_strt, indx_fut_end, dummy = _fetch_wthrset_indices(climgen.fut_wthr_set_defn,
                                                                            next_strt_yr, sim_end_yr)

    pettmp_sim = {}
    for metric in pettmp_hist:
        if indx_hist_end is not None:
            if indx_hist_end == -1:
                hist_seg = pettmp_hist[metric][indx_hist_strt:]
            else:
                hist_seg = pettmp_hist[metric][indx_hist_strt:indx_hist_end]

            pettmp_sim[metric] = hist_seg
            del hist_seg

        if indx_fut_end is not None:
            if indx_fut_end == -1:
                fut_seg = pettmp_fut[metric][indx_fut_strt:]
            else:
                fut_seg = pettmp_fut[metric][indx_fut_strt:indx_fut_end]

            pettmp_sim[metric] += fut_seg
            del fut_seg

    pettmp_sim = _apply_start_year_correction(sim_strt_yr, climgen.hist_wthr_set_defn, pettmp_sim)

    return pettmp_sim

def open_wthr_NC_sets(climgen):
    """
    C
    """
    hist_wthr_dsets = {}
    fut_wthr_dsets = {}

    for metric, ds_fname in zip(list(['precip', 'tas']), list(['ds_precip', 'ds_tas'])):
        hist_wthr_dsets[metric] = Dataset(climgen.hist_wthr_set_defn[ds_fname])
        fut_wthr_dsets[metric] = Dataset(climgen.fut_wthr_set_defn[ds_fname])

    return hist_wthr_dsets, fut_wthr_dsets

def fetch_WrldClim_data(lgr, lat, lon, climgen, nc_dsets, lat_indx, lon_indx, hist_flag=False):
    """
    C
    """
    filterwarnings("error")

    pettmp = {}
    for metric in list(['precip', 'tas']):
        if hist_flag:
            varname = climgen.hist_wthr_set_defn[metric]
            try:
                slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]
            except BaseException as err:
                print(ERROR_STR + str(err))
                return None
        else:
            varname = climgen.fut_wthr_set_defn[metric]
            slice = nc_dsets[metric].variables[varname][:, lat_indx, lon_indx]

        # test to see if cell data is valid, if not then this location is probably sea
        # =============================================================================
        if type(slice[0]) is MaskedConstant:
            pettmp = None
            mess = 'No data at lat: {} {}\tlon: {} {}\thist_flag: {}\n'.format(lat, lat_indx, lon, lon_indx, hist_flag)
            lgr.info(mess)
            # print(mess)
        else:
            pettmp[metric] = [float(val) for val in slice]

    return pettmp

def get_wthr_nc_coords(dset_defn, latitude, longitude):
    """
    C
    """
    lon_frst = dset_defn['lon_frst']
    lat_frst = dset_defn['lat_frst']
    resol_lat = dset_defn['resol_lat']
    resol_lon = dset_defn['resol_lon']
    max_lat_indx = len(dset_defn['latitudes']) - 1
    max_lon_indx = len(dset_defn['longitudes']) - 1

    lat_indx = int(round((latitude - lat_frst)/resol_lat))
    lon_indx = int(round((longitude - lon_frst)/resol_lon))

    if lat_indx < 0 or lat_indx > max_lat_indx:
        print(WARNING + 'latitude index {} out of bounds for latitude {}\tmax indx: {}'
                                                            .format(lat_indx, round(latitude, 4), max_lat_indx))
        return -1, -1

    if lon_indx < 0 or lon_indx > max_lon_indx:
        print(WARNING + 'longitude index {} out of bounds for longitude {}\tmax indx: {}'
                                                            .format(lon_indx, round(longitude, 4), max_lon_indx))
        return -1, -1

    return lat_indx, lon_indx

def check_clim_nc_limits(form, bbox_aoi = None, wthr_rsrce = 'CRU') :
    """
    this function checks that the specified bounding box lies within extent of the requested weather dataset
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True
    '''
    if hasattr(form, 'combo10w'):
        wthr_rsrce = form.combo10w.currentText()
    
    if wthr_rsrce == 'NASA' or wthr_rsrce == 'CRU':
        return limits_ok_flag
    '''
    lon_ll_aoi = float(form.w_ll_lon.text())
    lat_ll_aoi = float(form.w_ll_lat.text())
    lon_ur_aoi = float(form.w_ur_lon.text())
    lat_ur_aoi = float(form.w_ur_lat.text())

    wthr_rsrce = wthr_rsrce + '_hist'      # was + '_Day'
    lat_ur_dset = form.wthr_sets[wthr_rsrce]['lat_ur']
    lon_ur_dset = form.wthr_sets[wthr_rsrce]['lon_ur']
    lat_ll_dset = form.wthr_sets[wthr_rsrce]['lat_ll']
    lon_ll_dset = form.wthr_sets[wthr_rsrce]['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + wthr_rsrce + ' weather datasets')
    else:
        print('AOI lies outwith ' + wthr_rsrce + ' weather datasets - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag

def update_fetch_progress(last_time, nmasked, ncompleted, ncells):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:
        percnt_nremain = round(100 * ((ncells - ncompleted)/ncells), 2)

        scmplt = format(ncompleted, ',')
        smask = format(nmasked, ',')
        mess = '\rCells:  Completed: ' + scmplt
        mess += '\tMasked: ' + smask
        mess += '\tPercent: ' + str(percnt_nremain)
        stdout.flush()
        stdout.write(mess)
        last_time = new_time

    return last_time

def associate_climate(site_rec, climgen, pettmp_hist, pettmp_fut, report_flag=True):
    """
    this function associates each soil grid point with the most proximate climate data grid cell
    at the time of writing (Dec 2015) HWSD soil data is on a 30 arc second grid whereas climate data is on 30 or 15 or
     7.5 arc minute grid i.e. 0.5 or 0.25 or 0.125 of a degree
    """
    func_name =  __prog__ + ' associate_climate'

    proximate_keys = {}
    gran_lat_cell, gran_lon_cell, latitude, longitude = site_rec[:4]
    metric_list = pettmp_fut.keys()

    # TODO: find a more elegant methodology
    # =====================================
    for lookup_key in pettmp_hist['precip']:

        if lookup_key in pettmp_fut['precip']:
            if pettmp_fut['precip'][lookup_key] == None or pettmp_fut['tas'][lookup_key] == None or \
                        pettmp_hist['precip'][lookup_key] == None or pettmp_fut['tas'][lookup_key] == None:
                continue
            else:
                slat, slon = lookup_key.split('_')
                gran_lat = int(slat)
                gran_lon = int(slon)

                # situation where grid cell is coincidental with weather cell
                # ===========================================================
                if gran_lat == gran_lat_cell and gran_lon == gran_lon_cell:
                    climgen.lgr.info('Cell with lookup key ' + lookup_key + ' is coincidental with weather cell')
                    pettmp_out = {}
                    for metric in pettmp_fut.keys():
                        pettmp_out[metric] = list([pettmp_hist[metric][lookup_key], pettmp_fut[metric][lookup_key]])
                    return pettmp_out
                else:
                    proximate_keys[lookup_key] = list([gran_lat, gran_lon])

    # return empty dict if no proximate keys (unlikely)
    # =================================================
    if len(proximate_keys) == 0:
        if report_flag:
            mess = '\nNo weather keys assigned for site record with granular coordinates: '
            mess += '{} {}\tand lat/lon:'.format(round(gran_lat_cell, 4), round(gran_lon_cell, 4))
            print(mess + '{} {}'.format(round(latitude,4), round(longitude,4)))
        return {}

    # use the minimum distance to assign weather for specified grid cell
    # first stanza: calculate the squares of the distances in granular units between the grid cell and weathers cells
    # ===============================================================================================================
    dist = {}
    total_dist = 0
    for lookup_key in proximate_keys:
        gran_lat, gran_lon = proximate_keys[lookup_key]
        # _write_coords_for_key('\t', proximate_keys, lookup_key, func_name)

        dist[lookup_key] = (gran_lat - gran_lat_cell)**2 + (gran_lon - gran_lon_cell)**2
        total_dist += dist[lookup_key]

    # find key corresponding to the minimum value using conversions to lists
    # ======================================================================
    minval = sorted(dist.values())[0]
    lookup_key = list(dist.keys())[list(dist.values()).index(minval)]
    _write_coords_for_key('Selected weather key', climgen, proximate_keys, lookup_key, func_name)

    # construct result
    # ================
    pettmp_new_hist = {}
    pettmp_new_fut = {}
    for metric in metric_list:
        pettmp_new_hist[metric] = pettmp_hist[metric][lookup_key]
        pettmp_new_fut[metric] = pettmp_fut[metric][lookup_key]

    return (pettmp_new_hist, pettmp_new_fut)

def _write_coords_for_key(mess, climgen, proximate_keys, lookup_key, func_name):
    """
     C
    """
    # should go in log file - TODO
    gran_lat, gran_lon = proximate_keys[lookup_key]
    latitude  = 90 - gran_lat/GRANULARITY
    longitude = gran_lon/GRANULARITY - 180
    mess += ' Lat: {}\tGran lat: {}\tLon: {}\tGran lon: {}\tin {}'.format(latitude, gran_lat, longitude, gran_lon,
                                                                                                        func_name)
    climgen.lgr.info(mess)

    return
