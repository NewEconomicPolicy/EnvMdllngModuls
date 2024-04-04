"""
#-------------------------------------------------------------------------------
# Name:        prepareEcosseFiles.py
# Purpose:
# Author:      s03mm5
# Created:     08/12/2015
# Copyright:   (c) s03mm5 2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""
__version__ = '1.0.00'
__prog__ = 'prepare_ecosse_files.py'

# Version history
# ---------------
#
from os.path import normpath, join, lexists, basename
from os import makedirs
import csv
from time import time, sleep
from sys import stdout
from shutil import copyfile
from PyQt5.QtWidgets import QApplication

from thornthwaite import thornthwaite
from weather_datasets import write_csv_wthr_file
from glbl_ecss_cmmn_funcs import write_kml_file, write_manifest_file, input_txt_line_layout, write_signature_file

sleepTime = 5
GRANULARITY = 120
WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '

def _wthr_for_simulation(wthr_sets, climgen, pettmp_hist, pettmp_fut):
    """
    return spliced weather for simulation period
    """
    sim_start_year = climgen.sim_start_year
    sim_end_year = climgen.sim_end_year

    # must improve on this TODO
    # =========================
    wthr_rsrc = climgen.wthr_rsrc
    if wthr_rsrc == 'CHESS':
        hist_start_year = wthr_sets['CHESS_historic']['year_start']
        hist_end_year = wthr_sets['CHESS_historic']['year_end']
        fut_start_year = wthr_sets['chess_rcp26_01']['year_start']
    else:
        hist_start_year = wthr_sets['CRU_hist']['year_start']
        hist_end_year = wthr_sets['CRU_hist']['year_end']
        fut_start_year = wthr_sets['ClimGen_A1B']['year_start']

    pettmp_sim = {}
    if sim_start_year >= fut_start_year:
        indx_strt = 12*(sim_start_year - fut_start_year)
        for metric in pettmp_fut:
            pettmp_sim[metric] = pettmp_fut[metric][indx_strt:]
    else:
        # historic takes priority over future weather
        # ===========================================
        indx_hist_strt = 12*(sim_start_year - hist_start_year)
        indx_fut_strt  = 12*(hist_end_year + 1 - fut_start_year)
        for metric in pettmp_fut:
            pettmp_sim[metric] = pettmp_hist[metric][indx_hist_strt:] + pettmp_fut[metric][indx_fut_strt:]

    return pettmp_sim

def _make_met_files(clim_dir, latitude, climgen, pettmp_grid_cell):
    """
    feed annual temperatures to Thornthwaite equations to estimate Potential Evapotranspiration [mm/month]
    """
    func_name = __prog__ +  '_make_met_files'

    if not lexists(clim_dir):
        try:
            makedirs(clim_dir)
        except FileNotFoundError as err:
            print(ERROR_STR + str(err) + '\n\tin ' + func_name)
            sleep(sleepTime)
            exit(0)

    location = ' latitude: {}\tclimate directory: {}'.format(latitude, clim_dir)

    # ==============
    start_year = climgen.sim_start_year
    end_year   = climgen.sim_end_year
    precip = pettmp_grid_cell['precipitation']         #
    temper = pettmp_grid_cell['temperature']

    indx1 = 0
    for year in range(start_year, end_year + 1):
        fname = 'met{0}s.txt'.format(year)
        met_path = join(clim_dir, fname)

        indx2 = indx1 + 12

        precipitation = precip[indx1:indx2]            #
        tmean = temper[indx1:indx2]
        try:
            tmp_max = max(tmean)
        except ValueError as err:
            mess = WARN_STR + str(err) + ' for ' + location
            print(mess)
            break

        # PET
        # ===
        if tmp_max > 0.0:
            pet = thornthwaite(tmean, latitude, year)
        else:
            pet = [0.0]*12
            mess = WARN_STR + 'monthly temperatures all below zero for ' + location
            print(mess)

        # TODO: do something about occasional runtime warning...
        pot_evapotrans = [round(p, 2) for p in pet]
        precip_out     = [round(p, 2) for p in precipitation]
        tmean_out      = [round(t, 2) for t in tmean]

        # write file
        # ==========
        output = []
        for tstep, mean_temp in enumerate(tmean_out):
            output.append([tstep+1, precip_out[tstep], pot_evapotrans[tstep], mean_temp])

        with open(met_path, 'w', newline='') as fpout:
            writer = csv.writer(fpout, delimiter='\t')
            writer.writerows(output)

        indx1 += 12

    return

def make_ecosse_file(form, climgen, ltd_data, site_rec, province, pettmp_grid_cell):
    """
    generate sets of Ecosse files for each site
    where each site has one or more soils and each soil can have one or more dominant soils
    pettmp_grid_cell is climate data for this soil grid point
    """
    func_name = 'make_ecosse_file'

    gran_lat, gran_lon, latitude, longitude, area, mu_globals_props = site_rec
    sims_dir = form.sims_dir
    fut_clim_scen = climgen.fut_clim_scen

    # unpack weather
    # ==============
    pettmp_hist = {}
    pettmp_fut  = {}
    for metric in list(['precipitation','temperature']):
        pettmp_hist[metric] = pettmp_grid_cell[metric][0]
        pettmp_fut[metric]  = pettmp_grid_cell[metric][1]

    # calculate historic average weather
    # ==================================
    wthr_rsrc = climgen.wthr_rsrc
    if wthr_rsrc == 'CHESS':
        dset_start_year = form.wthr_sets['CHESS_historic']['year_start']
    else:
        dset_start_year = form.wthr_sets['CRU_hist']['year_start']

    hist_start_year = climgen.hist_start_year
    indx_start = 12*(hist_start_year - dset_start_year)

    hist_end_year = climgen.hist_end_year
    indx_end   = 12*(hist_end_year - dset_start_year + 1) # end year includes all 12 months - TODO: check

    # use dict-comprehension to initialise precip. and temperature dictionaries
    # =========================================================================
    hist_precip = {mnth: 0.0 for mnth in climgen.months}
    hist_tmean  = {mnth: 0.0 for mnth in climgen.months}

    for indx in range(indx_start, indx_end, 12):

        for imnth, month in enumerate(climgen.months):
            try:
                hist_precip[month] += pettmp_hist['precipitation'][indx + imnth]
                hist_tmean[month]  += pettmp_hist['temperature'][indx + imnth]
            except IndexError as err:
                mess = str(err) + ' in ' + func_name
                print(mess); form.lgr.info(mess)
                return

    # write stanza for input.txt file consisting of long term average climate
    # =======================================================================
    hist_wthr_recs = []
    num_hist_years = hist_end_year - hist_start_year + 1
    for month in climgen.months:
        ave_precip = hist_precip[month]/num_hist_years
        hist_wthr_recs.append(input_txt_line_layout('{}'.format(round(ave_precip,1)), \
                                            '{} long term average monthly precipitation [mm]'.format(month)))

    for month in climgen.months:
        ave_tmean = hist_tmean[month]/num_hist_years
        hist_wthr_recs.append(input_txt_line_layout('{}'.format(round(ave_tmean,2)), \
                                            '{} long term average monthly temperature [degC]'.format(month)))

    # write met files set for all simulations for this grid cell
    # ==========================================================
    study = form.study
    gran_coord = '{0:0=5g}_{1:0=5g}'.format(gran_lat, gran_lon)
    met_rel_path = '..\\' + gran_coord + '\\'
    clim_dir = normpath(join(sims_dir, study, gran_coord))
    simulation_wthr = _wthr_for_simulation(form.wthr_sets, climgen, pettmp_hist, pettmp_fut)
    _make_met_files(clim_dir, latitude, climgen, simulation_wthr)

    # create additional weather related files from already existing met files
    # =======================================================================
    irc = climgen.create_FutureAverages(clim_dir, latitude, gran_lat, longitude, gran_lon)

    # TODO: improve
    write_csv_wthr_file(form.lgr, study, wthr_rsrc, fut_clim_scen, latitude, longitude,
                                climgen.sim_start_year, climgen.sim_end_year,
                                pettmp_fut['precipitation'], pettmp_fut['temperature'], clim_dir)

    #------------------------------------------------------------------
    # Create a set of simulation input files for each dominant
    # soil-land use type combination
    #------------------------------------------------------------------
    # construct directory name with all dominant soils

    for pair in mu_globals_props.items():
        mu_global, proportion = pair
        area_for_soil = area*proportion
        soil_list = form.hwsd_mu_globals.soil_recs[mu_global]

        for soil_num, soil in enumerate(soil_list):
            identifer = 'lat{0:0=7d}_lon{1:0=7d}_mu{2:0=5d}_s{3:0=2d}'.format(gran_lat, gran_lon, mu_global, soil_num + 1)

            sim_dir = join(sims_dir, study, identifer)

            if not lexists(sim_dir):
                makedirs(sim_dir)

            ltd_data.write(sim_dir, soil, latitude, hist_wthr_recs, met_rel_path)

            # write kml file if requested and signature file
            # ==============================================
            if form.sttngs['kml_flag'] and soil_num == 0:
                write_kml_file(sim_dir,  str(mu_global), mu_global, latitude, longitude)

            write_signature_file(sim_dir, mu_global, soil, latitude, longitude, province)

            # copy across Model_Switches.dat file
            # ===================================
            outMdlSwtchs = join(sim_dir, basename(form.default_model_switches))
            copyfile(form.default_model_switches, outMdlSwtchs)

        # manifest file is essential for subsequent processing
        # ====================================================
        write_manifest_file(form.study, fut_clim_scen, sim_dir, soil_list, mu_global, latitude, longitude, area_for_soil)

    # end of Soil loop
    # ================

    return

def update_progress(last_time, completed, est_num_sims, skipped, warning_count, w_prgrss = None):
    """
    Update progress bar
    """
    new_time = time()
    if new_time - last_time > 5:

        if w_prgrss is None:
            mess = '\rCompleted: {:=6d} Skipped: {:=5d} Warnings: {:=5d} Remaining: {:=6d}' \
                .format(completed, skipped, warning_count, est_num_sims - completed)
            stdout.flush()
            stdout.write(mess)
        else:
            mess = '\rprocessed {:=6d} cells of which {:=6d} have been skipped with {:=6d} % remaining' \
                .format(completed,  skipped, est_num_sims - completed)
            w_prgrss.setText(mess)

        last_time = new_time
        QApplication.processEvents()

    return last_time
