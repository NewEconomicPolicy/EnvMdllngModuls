#-------------------------------------------------------------------------------
# Name:        getClimGenFns.py
# Purpose:     additional functions for getClimGenNC.py
# Author:      s03mm5
# Created:     08/02/2018
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'getClimGenFns.py'
__author__ = 's03mm5'

from time import time
from sys import stdout
from PyQt5.QtWidgets import QApplication

METRICS = ['precip', 'tas']

def associate_climate(site_rec, climgen, pettmp_hist, pettmp_fut):
    """
    this function associates each soil grid point with the most proximate climate data grid cell
    at the time of writing (Dec 2015) HWSD soil data is on a 30 arc second grid whereas climate data is on 30 or 15 or
     7.5 arc minute grid i.e. 0.5 or 0.25 or 0.125 of a degree
    """
    func_name =  __prog__ + ' associate_climate'

    proximate_keys = {}
    gran_lat_cell, gran_lon_cell, latitude, longitude, dummy, dummy = site_rec
    metric_list = pettmp_fut.keys()

    # TODO: find a more elegant methodology
    # =====================================
    for lookup_key in pettmp_hist['precipitation']:
        if pettmp_fut['precipitation'][lookup_key] == None or pettmp_fut['temperature'][lookup_key] == None or \
                    pettmp_hist['precipitation'][lookup_key] == None or pettmp_fut['temperature'][lookup_key] == None:
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
        print('\nNo weather keys assigned for site record with granular coordinates: {} {}\tand lat/lon: {} {}'
                                        .format(gran_lat_cell, gran_lon_cell, round(latitude,4), round(longitude,4)))
        return {}

    '''
    use the minimum distance to assign weather for specified grid cell
    '''

    # first stanza: calculate the squares of the distances in granular units between the grid cell and weathers cells
    # =============
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
    # _write_coords_for_key('Selected weather key', climgen, proximate_keys, lookup_key, func_name)

    #
    pettmp_final = {}
    for metric in metric_list:
        pettmp_final[metric] = list([pettmp_hist[metric][lookup_key], pettmp_fut[metric][lookup_key]])

    return pettmp_final

def check_clim_nc_limits(form, wthr_rsrc, bbox_aoi = None):

    """
    this function makes sure that the specified bounding box lies within extent of the requested weather dataset
    NB lats run from North to South
        lons run from West to East
    """
    func_name =  __prog__ + ' check_clim_nc_limits'

    limits_ok_flag = True

    # CRU is global
    # =============
    if wthr_rsrc == 'CRU':
        return limits_ok_flag
    #
    wthr_set_name = form.wthr_set_linkages[wthr_rsrc][0]
    wthr_set = form.wthr_sets[wthr_set_name]

    lon_ur_aoi = float(form.w_ur_lon.text())
    lat_ur_aoi = float(form.w_ur_lat.text())
    if form.version == 'HWSD_grid':
        lon_ll_aoi = float(form.w_ll_lon.text())
        lat_ll_aoi = float(form.w_ll_lat.text())
    else:
        lon_ll_aoi = lon_ur_aoi
        lat_ll_aoi = lat_ur_aoi

    lat_ur_dset = wthr_set['lat_ur']
    lon_ur_dset = wthr_set['lon_ur']
    lat_ll_dset = wthr_set['lat_ll']
    lon_ll_dset = wthr_set['lon_ll']

    # similar functionality in lu_extract_fns.py in LU_extract project
    # ================================================================
    if (lon_ll_dset < lon_ll_aoi and lon_ur_dset > lon_ur_aoi) and \
                    (lat_ll_dset < lat_ll_aoi and lat_ur_dset > lat_ur_aoi):
        print('AOI lies within ' + wthr_rsrc + ' weather dataset')
    else:
        print('AOI lies outwith ' + wthr_rsrc + ' weather dataset - LL long/lat: {} {}\tUR long/lat: {} {}'
              .format(lon_ll_dset, lat_ll_dset, lon_ur_dset, lat_ur_dset))
        limits_ok_flag = False

    return limits_ok_flag

def update_progress_clim_soil(last_time, nsoilres, pt_key, ncsv_lines, skipped = 0, failed = 0):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:

        mess = '\rSize of soil list: {}\tpt_key: {}\tNumber of sites remaining: {}'\
                                            .format(nsoilres, pt_key, ncsv_lines - nsoilres)
        stdout.flush()
        stdout.write(mess)
        last_time = new_time

    return last_time

def update_progress_clim(metric, last_time, total_num, ngrid_cells, no_data, w_prgrss = None):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 2:

        if w_prgrss is None:
            mess = '\rCompleted checking of: {:=7d} climate cells\tNo data: {:=7d}\tRemaining: {:=7d}\tTotal: {:=7d}' \
                .format(ngrid_cells, no_data, 1 + total_num - ngrid_cells, total_num)
            stdout.flush()
            stdout.write(mess)
        else:
            mess = '\rChecked {:=7d} {:6s} climate cells of which {:=7d} have no data   {:=5.2f} % complete' \
                .format(ngrid_cells, metric, no_data, 100*ngrid_cells/total_num)
            w_prgrss.setText(mess)

        last_time = new_time
        QApplication.processEvents()

    return last_time
