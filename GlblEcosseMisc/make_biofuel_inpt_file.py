#-------------------------------------------------------------------------------
# Name:
# Purpose:     make and input file for biofuel program from CHESS data
# Author:      Mike Martin
# Created:     18 June 2023
# Description:
#-------------------------------------------------------------------------------
#
__author__ = 's03mm5'
__prog__ = 'make_biofuel_inpt_file.py'
__version__ = '0.0'

from os import mkdir, getcwd
from os.path import isfile, isdir, join, normpath, split
from sys import exit, stdout
from time import sleep, time
from pandas import read_excel, read_csv
from csv import writer

from weather_datasets_lite import read_wthr_dsets_detail
from dset_fns_and_class import ElevationSet, open_proj_NC_sets, close_proj_NC_sets, fetch_wthr_elev

sleepTime = 5
PROGRAM_ID = 'biofuel_inpt_file'
ERROR_STR = '*** Error *** '

RQURD_WTHR_RSRCS = 'CHESS'
SCENARIOS = list(['rcp26', 'rcp45', 'rcp60', 'rcp85'])
REALISATIONS = list(['01', '04', '06', '15'])

MAX_CELLS = 200     # limit on output

class Form(object):
    """
    Class to permit bulk monthly file generation in batch mode
    """
    def __init__(self, datasets_dir):
        """
        make monthly directories is required
        """
        wthr_dir = datasets_dir
        if not isdir(wthr_dir):
            print(ERROR_STR + 'weather directory ' + wthr_dir + ' must exist')
            sleep(sleepTime)
            exit(0)

        elev_fn = join(datasets_dir, 'elevation', 'elevation_UK.nc')
        if not isfile(elev_fn):
            print(ERROR_STR + 'elevation dataset ' + elev_fn + ' must exist')
            sleep(sleepTime)
            exit(0)

        '''
        grid_ref_xls = join(datasets_dir, 'grid_refs', 'UK_grid_refs_1km.xlsx')
        if not isfile(grid_ref_xls):
            print(ERROR_STR + 'grid references Excel file ' + grid_ref_xls + ' must exist')
            sleep(sleepTime)
            exit(0)        
        '''
        grid_ref_csv = join(datasets_dir, 'grid_refs', 'UK_grid_refs_1km.csv')
        if not isfile(grid_ref_csv):
            print(ERROR_STR + 'grid references CSV file ' + grid_ref_csv + ' must exist')
            sleep(sleepTime)
            exit(0)

        rslt_dir = join(datasets_dir, 'result')
        if not isdir(rslt_dir):
            mkdir(rslt_dir)

        self.sttngs = {'wthr_dir': wthr_dir, 'elev_fn': elev_fn, 'grid_ref_csv': grid_ref_csv, 'rslt_dir': rslt_dir}

        cw_dir = getcwd()
        self.fobj_masked = open(cw_dir + '\\masked_cells.txt', 'w')

        self.fobj_inpts = open(cw_dir + '\\biofuel_inputs.txt', 'w', newline='')
        self.writer = writer(self.output_fhs, delimiter=',')

        hdr_rec = 'junk'
        self.writer.writerow(hdr_rec)

    def make_input_file(self):
        """
         
        """
        resol = 10 * 1000   # metres

        read_wthr_dsets_detail(self, RQURD_WTHR_RSRCS)
        wthr_defn = self.wthr_sets['CHESS_historic']

        grid_ref_csv = self.sttngs['grid_ref_csv']
        grid_ref_df = read_csv(grid_ref_csv, sep = ',')

        elev_fn = self.sttngs['elev_fn']
        elev_defn = ElevationSet(elev_fn)

        open_proj_NC_sets(elev_defn, wthr_defn)

        ic = 0
        for estng, nrthng in zip(grid_ref_df['Grid_Easting'], grid_ref_df['Grid_Northing']):
            if (estng - 500) % resol == 0 and (nrthng - 500) % resol == 0:
                pettmp, elev, lat, lon = fetch_wthr_elev(self.fobj_masked, estng, nrthng, elev_defn, wthr_defn)
                if pettmp is None:
                    continue

                print('Data at E/N: {} {}'.format(estng, nrthng))
                self.writer.writerow('record')


        close_proj_NC_sets(elev_defn, wthr_defn)
        self.fobj.close()

        print('Finished')

        return True

def _get_year(nc_fname):
    '''

    '''
    short_daily_fn = split(nc_fname)[1]
    fn_cmpnts = short_daily_fn.split('monthly')
    this_year = fn_cmpnts[1][1:5]

    return this_year

def main_biofuel():
    """
    Entry point
    """
    datasets_dir = normpath('E:\\Faith_Sadiq\\biofuel_datasets')

    form = Form(datasets_dir)

    form.make_input_file()

    print('\nfinished conversion')     # cosmetic

    return
