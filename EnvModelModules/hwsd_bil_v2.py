#-------------------------------------------------------------------------------
# Name:        hwsd_bil_v2.py
# Purpose:     Class to read HWSD .hdr and .bil files and return a 2D array comprising MU_GLOBAL values for a given
#              bounding box
#               input files are: hwsd.hdr and hwsd.bil
#               notes:
#                   a) lat longs are stored in unites of 30 arc seconds for sake of accuracy and simplicity - so we
#                       convert lat/longs read from shapefiles to nearest seconds using round()
#
#                   b) output file names take the form: mu_global_NNNNN.nc where NNNNN is the region name??
#
#                   c) modified on 03-09-2016 to include soils which have only a top layer
#
# Author:      Mike Martin
# Created:     16/09/2025
# Licence:     <your licence>
#-------------------------------------------------------------------------------

__prog__ = 'hwsd_bil_v2.py'
__version__ = '0.0.0'

from os.path import join, exists, isdir,isfile
from pandas import read_sql
from sys import exit, stdout
from struct import unpack
from locale import LC_ALL, setlocale, format_string
from numpy import arange, dtype, zeros, int32
from time import sleep
from pyodbc import connect, drivers
from tabulate import tabulate

ERROR_STR = '*** Error *** '

VAR_FLOAT_LIST = ['ulxmap', 'ulymap', 'xdim', 'ydim']
VAR_STR_LIST = ['pixeltype', 'byteorder', 'layout']

sleepTime = 5

def fetch_metadata(cursor):
    """
    C
    """
    cmd = 'select * from HWSD2_LAYERS_METADATA'
    try:
        # df = read_sql(cmd, conn)
        cursor.execute(cmd)
    except BaseException as err:
        print(str(err))
        return

    recs_lyrs = [rec for rec in cursor.fetchall()]

    cmd = 'select * from HWSD2_SMU_METADATA'
    cursor.execute(cmd)
    recs_smu = [rec for rec in cursor.fetchall()]

    return recs_lyrs, recs_smu

def fetch_accesss_conn(hwsd_dir):
    """
    C
    """
    ms_drvr = 'Microsoft Access Driver(*.mdb, *.accdb)'
    access_db_fn = join(hwsd_dir, 'mdb', 'HWSD2.mdb')

    ms_srch_str = 'Microsoft Access Driver'
    drvr_nms = [drvr_nm for drvr_nm in drivers() if drvr_nm.startswith(ms_srch_str)]
    if len(drvr_nms) == 0:
        print(ERROR_STR + 'could not find ' + ms_srch_str + ' among ODBC drivers')
        return

    ms_drvr = drvr_nms[0]

    conn = connect(Driver=ms_drvr, DBQ=access_db_fn)
    return conn, conn.cursor()

def _make_four_line_table(coverage, mu_global, wrb2_value, wrb2, fao90_value, fao90):
    """
    make first four lines
    """

    a = [
        ['Coverage:', coverage],
        ['Soil Mapping Unit (SMU):', mu_global],
        ['Dominant Soil Unit (WRB 2022):', wrb2_value + '(' + wrb2 + ')'],
        ['Dominant Soil Unit (FAO 1990):', fao90_value + '(' + fao90 + ')']
    ]
    # create header
    headers = ['Name', 'City']

    print(tabulate(a, headers=headers, tablefmt='grid'))

    return

def get_soil_recs(cursor, mu_globals):
    """
    # retcode = cursor.execute('select * from D_ROOTS')
    """
    table_names = [table_info.table_name for table_info in cursor.tables(tableType='TABLE')]
    layer_column_names = [row.column_name for row in cursor.columns(table='HWSD2_LAYERS')]
    mu_globals = [mu_global for mu_global in mu_globals.keys()]
    mu_global = mu_globals[0]

    # fetch Coverage, Dominant Soil Units for WRB and FAO
    # ===================================================
    VARS = ' COVERAGE, WRB2, FAO90 '
    cmd = 'select ' + VARS + ' from HWSD2_SMU where HWSD2_SMU_ID = ' + str(mu_global)
    cursor.execute(cmd)

    smu_recs = [rec for rec in cursor.fetchall()]
    code, wrb2, fao90 = smu_recs[0]

    cmd = 'select VALUE from D_COVERAGE where CODE = ' + str(code)
    cursor.execute(cmd)
    recs = [rec for rec in cursor.fetchall()]
    coverage = recs[0][0]

    cmd = "select VALUE from D_WRB2 where CODE = '" + wrb2 + "'"
    cursor.execute(cmd)

    recs = [rec for rec in cursor.fetchall()]
    wrb2_value = recs[0][0]

    cmd = "select VALUE from D_FAO90 where CODE = '" + fao90 + "'"
    cursor.execute(cmd)
    recs = [rec for rec in cursor.fetchall()]
    fao90_value = recs[0][0]

    _make_four_line_table(coverage, mu_global, wrb2_value, wrb2, fao90_value, fao90)

    VARS = ' SEQUENCE, SHARE, LAYER, SAND, SILT, CLAY, BULK, REF_BULK, ORG_CARBON, PH_WATER '
    cmd = 'select ' + VARS + ' from HWSD2_LAYERS where HWSD2_SMU_ID = ' + str(mu_global)

    # layer_df = read_sql(cmd, conn)

    retcode = cursor.execute(cmd)
    layer_recs = [rec for rec in cursor.fetchall()]

    return layer_recs

def check_hwsd_integrity(hwsd_dir):
    """
    invoked at startup - check essential directory and files
    """
    integrity_flag = True

    if not isdir(hwsd_dir):
        print(ERROR_STR + 'HWSD directory ' + hwsd_dir + ' must exist')
        integrity_flag = False

    # check main table and header file
    # =================================
    inp_fname = 'HWSD2.mdb'
    mdb_file = join(hwsd_dir, 'mdb', inp_fname)
    if not isfile(mdb_file):
        print(ERROR_STR + 'HWSD database file HWSD table ' + mdb_file + ' must exist')
        integrity_flag = False

    inp_fname = 'HWSD2.hdr'
    hdr_file = join(hwsd_dir, 'raster', inp_fname)
    if not isfile(hdr_file):
        print(ERROR_STR + 'HWSD header file ' + hdr_file + ' must exist')
        integrity_flag = False

    if integrity_flag:
        print('HWSD version 2 in ' + hwsd_dir + ' has passed integrity check')
        return
    else:
        sleep(sleepTime)
        exit(0)

def validate_hwsd_rec (lgr, mu_global, data_rec):
    """
    # function to make sure we are supplying a valid HWSD data record for a given mu_global

    # TODO: improve this function by using a named tuple:
                 data_rec = list([
                 0 row['SHARE'],
                 1 row['T_SAND'],
                 2 row['T_SILT'],
                 3 row['T_CLAY'],
                 4 row['T_BULK_DENSITY'],
                 5 row['T_OC'],
                 6 row['T_PH_H2O'],
                 7 row['S_SAND'],
                 8 row['S_SILT'],
                 9 row['S_CLAY'],
                10 row['S_BULK_DENSITY'],
                11 row['S_OC'],
                12 row['S_PH_H2O']
    """
    share, t_sand, t_silt, t_clay, t_bulk_density, t_oc, t_ph_h2o, \
        s_sand, s_silt, s_clay, s_bulk_density, s_oc, s_ph_h2o = data_rec

    # topsoil is mandatory
    ts_rec = data_rec[1:7]
    for val in ts_rec:
        if val == '':
            # log a message
            mess =  'Incomplete topsoil data - mu_global: ' + str(mu_global) + '\trecord: ' + ', '.join(data_rec)
            # print(mess)
            lgr.info(mess)
            return False

    '''
    MJM: borrowed from mksims.py:
    Convert top and sub soil Organic Carbon percentage weight to kgC/ha as follows:
                      / 100: % to proportion,
                     * 100000000: cm3 to hectares,
                      /1000: g to kg,
                     *30: cm to lyr depth
    '''
    # lgr.info('Soil data check for mu_global {}\tshare: {}%'.format(mu_global,share))
    t_bulk = float(t_bulk_density)  # T_BULK_DENSITY
    t_oc = float(t_oc)    # T_OC
    t_c = round(t_bulk * t_oc / 100.0 * 30 * 100000000 / 1000.0, 1)

    # sub soil
    ss_rec = data_rec[7:]
    for val in ss_rec:
        if val == '':
            # log a message
            lgr.info('Incomplete subsoil data - mu_global: ' + str(mu_global) + '\trecord: ' + ', '.join(data_rec))
            # modified share to from int to float
            return list([t_c, t_bulk, float(t_ph_h2o), int(t_clay), int(t_silt), int(t_sand), float(share)])

    s_bulk = float(s_bulk_density) # S_BULK_DENSITY
    s_oc = float(s_oc)   # S_OC
    s_c = round(s_bulk * s_oc / 100.0 * 70 * 100000000 / 1000.0, 1)

    # Two layers: C content, bulk, PH, clay, silt, sand
    try:
        # modified share to from int to float
        ret_list = list([t_c, t_bulk, float(t_ph_h2o), int(t_clay), int(t_silt), int(t_sand),
                 s_c, s_bulk, float(s_ph_h2o), int(s_clay), int(s_silt), int(s_sand), float(share)])
    except ValueError as e:
        ret_list = []
        print('problem {} with mu_global {}\tdata_rec: {}'.format(e, mu_global, data_rec))
        sleep(sleepTime)
        exit(0)

    return ret_list

class HWSD_bil(object,):

    def __init__(self, lgr, hwsd_dir):
        """
        open header file and read content consisting 14 lines
        """
        inp_file = join(hwsd_dir, 'raster', 'HWSD2.hdr')
        with open(inp_file, 'r') as finp:
            lines = finp.readlines()
        lines = [line.rstrip('\n') for line in lines]

        hdr_dict = {}
        for line in lines:
            var, val = line.split()
            var = var.lower()

            if var in VAR_FLOAT_LIST:
                val = float(val)
            elif var in VAR_STR_LIST:
                pass
            else:
                val = int(val)

            hdr_dict[var] = val

        # the HWSD grid covers the globe's land area with 30 arc-sec grid-cells
        # =====================================================================
        self.granularity = int(1.0 / hdr_dict['xdim'])

        self.hwsd_dir = hwsd_dir
        self.lgr = lgr
        self.nrows = hdr_dict['nrows']
        self.ncols = hdr_dict['ncols']
        self.wordsize = round(hdr_dict['nbits']/8)      # wordsize is in bytes

        # these attrubutes change according to the bounding box
        # =====================================================
        self.rows = arange(1)   # create ndarray
        self.mu_globals = []
        self.badCntr = 0
        self.zeroCntr = 0

        # these define the extent of the grid
        # ===================================
        self.nlats = 0; self.nlons = 0
        self.nrow1 = 0; self.nrow2 = 0
        self.ncol1 = 0; self.ncol2 = 0

        '''
        this will store dictionary of mu_global keys and corresponding soil data
        and will be refreshed each time user selects a new bounding box
        '''
        self.soil_recs = {}

    def read_bbox_mu_globals(self, bbox, snglPntFlag = False):
        """
        this function creates a grid of MU_GLOBAL values corresponding to a given
        bounding box defined by two lat/lon pairs

        the HWSD grid covers the globe's land area with 30 arc-sec grid-cells
        AOI is typically county sized e.g. Isle of Man
        """
        func_name =  __prog__ + ' read_bbox_mu_globals'
        self.lgr.info('Running programme ' + func_name)

        granularity = self.granularity
        ncols = self.ncols
        wordsize = self.wordsize

        if snglPntFlag:
            nlats = 1
            nlons = 1
            lower_left_coord = [bbox[1],bbox[0]]
            nrow1 = round((90.0 - lower_left_coord[0])*granularity)
            nrow2 = nrow1
            ncol1 = round((180.0 + lower_left_coord[1])*granularity)
            ncol2 = ncol1
            nlats = 1
            nlons = 1
        else:
            # these are lat/lon
            upper_left_coord = [bbox[3],bbox[0]]
            lower_right_coord = [bbox[1],bbox[2]]
            self.upper_left_coord = upper_left_coord
            self.lower_right_coord = lower_right_coord

            nlats = round((upper_left_coord[0] - lower_right_coord[0])*granularity)
            nlons = round((lower_right_coord[1] - upper_left_coord[1])*granularity)
            #
            # work out first and last rows
            nrow1 = round((90.0 - upper_left_coord[0])*granularity) + 1
            nrow2 = round((90.0 - lower_right_coord[0])*granularity)
            ncol1 = round((180.0 + upper_left_coord[1])*granularity) + 1
            ncol2 = round((180.0 + lower_right_coord[1])*granularity)

        chunksize = wordsize*nlons
        form = str(int(chunksize/2)) + 'H'  # format for unpack

        # TODO - construct a 2D array from extraction
        inpfname = 'HWSD2.bil'
        inp_file = join(self.hwsd_dir, 'raster', inpfname)
        finp = open(inp_file, 'rb')

        nrows_read = 0
        rows = arange(nlats*nlons)
        rows.shape = (nlats,nlons)

        # read subset of HWSD .bil file
        # =============================
        for nrow in range(nrow1, nrow2 + 1):
            offst = (ncols*nrow + ncol1)*wordsize

            finp.seek(offst)
            chunk = finp.read(chunksize)
            if chunk:
                rows[nrows_read:] = unpack(form, chunk)
                # build array of mu_globals
                nrows_read += 1
                # print(' row {0:d}\tform {1:s}'.format(nrow,form))
            else:
                break

        finp.close()

        self.nrow1 = nrow1
        self.nrow2 = nrow2
        self.ncol1 = ncol1
        self.ncol2 = ncol2

        self.rows = rows
        self.nlats = nlats
        self.nlons = nlons

        self.lgr.info('Exiting function {0} after reading {1} records from .bil file\n'.format(func_name,nrows_read))
        return nrows_read

    def mu_global_list(self):
        """
        build a list of mu_globals from grid
        """
        # reshape
        mu_globals = []
        nlats = self.nlats
        nlons = self.nlons

        # reshape tuple
        self.rows.shape = (nlats*nlons)
        for mu_global in self.rows:
            if mu_global not in mu_globals:
                mu_globals.append(mu_global)

        # revert shape of tuple
        self.rows.shape = (nlats,nlons)

        self.mu_globals = mu_globals

        return len(mu_globals)

    def get_mu_globals_dict(self):

        func_name =  __prog__ + ' get_mu_globals_dict'
        # function returns list of MU_GLOBALs in grid with number of occurrences of each MU_GLOBAL value
        # i.e. build dictionary of MU_GLOBALs from grid

        mu_globals = {}     # initialise dictionary

        nlats = self.nlats
        nlons = self.nlons

        for irow in range(0,nlats):

            val_rec = []
            for icol in range(0,nlons):

                mu_global = self.rows[irow,icol]

                if mu_global not in mu_globals.keys():
                    # add detail
                    mu_globals[mu_global] = 1
                else:
                    mu_globals[mu_global] += 1

        self.lgr.info('Func: {0}\tUnique number of mu_globals: {1}'.format(func_name, len(mu_globals)))
        for key in sorted(mu_globals.keys()):
            self.lgr.info('\t' + str(key) + ' :\t' + str(mu_globals[key]))

        # filter condition where there is only sea, i.e. mu_global == 0
        if len(mu_globals.keys()) == 1:
            for key in mu_globals.keys():
                if key == 0:
                    mu_globals = {}
                    print('Only one mu_global with value of zero - nothing to process')

        return mu_globals

def fetch_accesss_cursor(hwsd_dir):
    """
    C
    """
    ms_drvr = 'Microsoft Access Driver(*.mdb, *.accdb)'
    access_db_fn = join(hwsd_dir, 'mdb', 'HWSD2.mdb')

    ms_srch_str = 'Microsoft Access Driver'
    drvr_nms = [drvr_nm for drvr_nm in drivers() if drvr_nm.startswith(ms_srch_str)]
    if len(drvr_nms) == 0:
        print(ERROR_STR + 'could not find ' + ms_srch_str + ' among ODBC drivers')
        return

    ms_drvr = drvr_nms[0]

    conn = connect(Driver=ms_drvr, DBQ=access_db_fn)
    return conn.cursor()

