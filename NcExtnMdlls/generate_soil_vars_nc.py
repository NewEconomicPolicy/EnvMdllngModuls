"""
#-------------------------------------------------------------------------------
# Name:        hwsd_glblecsse_fns.py
# Purpose:     consist of high level functions invoked by main GUI
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
# Description:
#   comprises two functions:
#       def _write_to_soil_files(form, climgen,  mask_defn, num_band)
#       def generate_banded_sims(form)
#-------------------------------------------------------------------------------
#
"""
__prog__ = 'generate_soil_vars_grid.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

import time
from netCDF4 import Dataset
from PyQt5.QtWidgets import QApplication

from netcdf_funcs import create_soil_nc_dset, get_nc_coords
from glbl_ecsse_high_level_sp import simplify_soil_recs, _simplify_aoi

import hwsd_bil

from hwsd_mu_globals_fns import gen_grid_cells_for_band
from prepare_ecosse_files import update_progress
from glbl_ecss_cmmn_cmpntsGUI import calculate_grid_cell

SOIL_METRICS = dict({'C_content': 'kgC/ha', 'bulk_dens': 'g/cm3', 'pH': 'pH', 'clay': '%',
                     'silt': '%', 'sand ': '%'})

def _write_to_soil_nc(form, nc_dset, max_lon_indx, max_lat_indx, num_band):
    """
    Main loop for generating soil data outputs
    """
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)   # instantiate a soil grid
    calculate_grid_cell(form, hwsd.granularity)  # add requested grid resolution attributes to the form object
    bbox = form.bbox

    # create grid of mu_globals based on bounding box
    # ===============================================
    nvals_read = hwsd.read_bbox_hwsd_mu_globals(bbox, form.hwsd_mu_globals, form.req_resol_upscale)

    # retrieve dictionary consisting of mu_globals (keys) and number of occurences (values)
    # =====================================================================================
    print('\nRetrieving soil data for band ' + str(num_band))
    QApplication.processEvents()

    mu_globals = hwsd.get_mu_globals_dict()
    if mu_globals is None:
        print('No soil records for AOI: {}\n'.format(bbox))
        return

    mess = 'Retrieved {} values  of HWSD grid consisting of {} rows and {} columns: ' \
          '\n\tnumber of unique mu_globals: {}'.format(nvals_read, hwsd.nlats, hwsd.nlons, len(mu_globals))
    form.lgr.info(mess)

    # for each grid point in the band) return a quintuples of integer Lat/Lon 30 seconds coords, Lat/Lons and mu_global
    hwsd.bad_muglobals = form.hwsd_mu_globals.bad_mu_globals
    aoi_res, bbox = gen_grid_cells_for_band(hwsd, form.req_resol_upscale)
    if form.w_use_high_cover.isChecked():
        aoi_res = _simplify_aoi(aoi_res)

    lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi = bbox
    num_meta_cells = len(aoi_res)
    print('Band aoi LL lon/lat: {} {}\tUR lon/lat: {} {}\t# meta cells: {}'
                            .format(lon_ll_aoi, lat_ll_aoi, lon_ur_aoi, lat_ur_aoi, num_meta_cells))
    QApplication.processEvents()

    if num_meta_cells == 0:
        mess = 'No aoi_res recs therefore unable to create simulation files... \n'
        print(mess); form.lgr.info(mess)
        return

    mess = 'Generated {} Area of Interest grid cells for band {} '.format(num_meta_cells, num_band)
    form.lgr.info(mess); print(mess)

    print('Writing soil data for band {}...'.format(num_band))
    QApplication.processEvents()

    # open land use NC dataset
    # ========================
    last_time = time.time()
    start_time = time.time()
    completed = 0
    skipped = 0
    warn_count = 0

    # each soil can have one or more dominant soils
    # =======================================================================
    for site_rec in aoi_res:
        gran_lat, gran_lon, lat, long, area, mu_globals_props = site_rec
        lat_indx, lon_indx = get_nc_coords(form.bbox_nc, lat, long, max_lat_indx, max_lon_indx)
        if lat_indx == -1:
            continue

        for pair in mu_globals_props.items():
            mu_global, proportion = pair
            soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

            # soil_c, bulk_dens, ph = soil[:3]
            # ===============================
            nc_dset.variables['mu_global'][lat_indx, lon_indx] = mu_global
            for soil in soil_list:
                for metric, val in zip(SOIL_METRICS, soil[:3]):
                    nc_dset.variables[metric][lat_indx, lon_indx, 0] = val

                if len(soil) == 13:
                    for metric, val in zip(SOIL_METRICS, soil[6:9]):
                        nc_dset.variables[metric][lat_indx, lon_indx, 1] = val

        completed += 1
        last_time = update_progress(last_time, start_time, completed, num_meta_cells, skipped, warn_count)
        QApplication.processEvents()

    mess = '\nBand: {}\tskipped: {}\tcompleted: {}'.format(num_band, skipped, completed)
    print(mess); QApplication.processEvents()

    print('')   # spacer
    return completed

def make_soil_nc_outputs(form):
    '''
    called from GUI
    '''
    if hasattr(form, 'w_max_cells'):
        max_cells = int(form.w_max_cells.text())
    else:
        max_cells = 100000000

    if form.hwsd_mu_globals == None:
        print('Undetermined HWSD aoi - please select a valid HSWD csv file')
        return

    if form.w_use_dom_soil.isChecked():
        use_dom_soil_flag = True
    else:
        use_dom_soil_flag = False

    # make sure bounding box is correctly set
    # =======================================
    lon_ll = -11.5
    lat_ll = 34.5
    lon_ur = 35.0
    lat_ur = 72.0
    form.bbox =  list([lon_ll, lat_ll, lon_ur, lat_ur])

    # lat_ll_aoi is the floor i.e. least latitude, of the HWSD aoi which marks the end of the banding loop
    # ====================================================================================================
    lat_ll_aoi = form.hwsd_mu_globals.lat_ll_aoi
    lon_ll_aoi = form.hwsd_mu_globals.lon_ll_aoi
    lat_ur_aoi = form.hwsd_mu_globals.lat_ur_aoi
    lon_ur_aoi = form.hwsd_mu_globals.lon_ur_aoi
    bbox_aoi = list([lon_ll_aoi,lat_ll_aoi,lon_ur_aoi,lat_ur_aoi])

    # check overlap - study too far to west or east or too far south or north of AOI file
    # ===================================================================================
    if (lon_ur < lon_ll_aoi) or (lon_ll > lon_ur_aoi) or  (lat_ur < lat_ll_aoi) or (lat_ll > lat_ur_aoi):
        print('Error: Study bounding box and HWSD CSV file do not overlap - no simulations are possible')
        return

    # ============================ for each PFT end =====================================
    # print('Study bounding box and HWSD CSV file overlap')
    #        ============================================
    start_at_band = form.start_at_band
    print('Starting at band {}'.format(start_at_band))

    # extract required values from the HWSD database and simplify if requested
    # ========================================================================
    hwsd = hwsd_bil.HWSD_bil(form.lgr, form.hwsd_dir)

    # TODO: patch to be sorted
    # ========================
    mu_global_pairs = {}
    for mu_global in form.hwsd_mu_globals.mu_global_list:
        mu_global_pairs[mu_global] = None

    soil_recs = hwsd.get_soil_recs(mu_global_pairs)  # list is already sorted

    # TODO: patch to be sorted
    # ========================
    for mu_global in hwsd.bad_muglobals:
        del(soil_recs[mu_global])

    form.hwsd_mu_globals.soil_recs = simplify_soil_recs(soil_recs, use_dom_soil_flag)
    form.hwsd_mu_globals.bad_mu_globals = [0] + hwsd.bad_muglobals
    del(hwsd); del(soil_recs)

    # NC file is deleted if it already exists
    # ====================================='
    nc_fname = create_soil_nc_dset(form, lon_ll, lat_ll, lon_ur, lat_ur, SOIL_METRICS)
    if isinstance(nc_fname, int):
        return
    try:
        nc_dset = Dataset(nc_fname, 'a', format='NETCDF4')
    except TypeError as err:
        print('Unable to open output file. {}'.format(err))
        return
    max_lat_indx = nc_dset.variables['latitude'].shape[0] - 1
    max_lon_indx = nc_dset.variables['longitude'].shape[0] - 1

    # main banding loop
    # =================
    ncompleted = 0
    lat_step = 0.5
    nsteps = int((lat_ur-lat_ll)/lat_step) + 1
    for isec in range(nsteps):
        lat_ll_new = lat_ur - lat_step
        num_band = isec + 1
        '''
        if num_band > 3:       # TODO remove when no longer needed
            print('Exiting from processing after {} bands'.format(num_band - 1))
            break
        '''
        # if the latitude floor of the band has not reached the ceiling of the HWSD aoi then skip this band
        if lat_ll_new > form.hwsd_mu_globals.lat_ur_aoi or num_band < start_at_band:
            print('Skipping out of area band {} of {} with latitude extent of min: {}\tmax: {}\n'
              .format(num_band, nsteps, round(lat_ll_new,6), round(lat_ur, 6)))
        else:
            print('\nProcessing band {} of {} with latitude extent of min: {}\tmax: {}'
                                                    .format(num_band, nsteps, round(lat_ll_new,6), round(lat_ur, 6)))
            QApplication.processEvents()

            form.bbox = list([lon_ll, lat_ll_new, lon_ur, lat_ur])
            completed = _write_to_soil_nc(form, nc_dset, max_lon_indx, max_lat_indx, num_band)  # does actual work
            if isinstance(completed, int):
                ncompleted += completed
                if ncompleted >= max_cells:
                    print('Finished processing after generating {} cells'.format(ncompleted))
                    break

        # check to see if the last band is completed
        if lat_ll_aoi > lat_ll_new or num_band == nsteps:
            print('Finished processing after {} bands of latitude extents'.format(num_band))
            for ichan in range(len(form.fstudy)):
                form.fstudy[ichan].close()
            break

        lat_ur = lat_ll_new

    # close nc file
    # =============
    nc_dset.close()
    print('\nFinished soil metric writing')

    return

