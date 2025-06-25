"""
#-------------------------------------------------------------------------------
# Name:        glbl_ecss_cmmn_cmpntsGUI.py
# Purpose:     module comprising miscellaneous functions common to some Global Ecosse variations
# Author:      Mike Martin
# Created:     31/07/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
"""

__prog__ = 'glbl_ecss_cmmn_cmpntsGUI.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QLabel, QComboBox, QLineEdit, QPushButton

from shape_funcs import calculate_area, format_bbox

sleepTime = 5

WARN_STR = '*** Warning *** '
ERROR_STR = '*** Error *** '
DUMMY_STR = 'dummy'

GRANULARITY = 120

WDGT_SIZE_60 = 60
WDGT_SIZE_80 = 80

RESOLUTIONS = {120:'30"', 30:'2\'', 20:'3\'', 10:'6\'', 8:'7\' 30"', 6:'10\'', 4:'15\'', 3:'20\'', 2:'30\''}
REVERSE_RESOLS = {}
for key in RESOLUTIONS:
    REVERSE_RESOLS[RESOLUTIONS[key]] = key

def glblecss_bounding_box(form, grid, irow):
    """

    """
    irow += 1
    icol = 0
    # UR lon/lat
    # ==========
    irow += 1
    lbl02a = QLabel('Upper right longitude:')
    lbl02a.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl02a, irow, 0)

    w_ur_lon = QLineEdit()
    w_ur_lon.setFixedWidth(WDGT_SIZE_80)
    grid.addWidget(w_ur_lon, irow, 1)
    form.w_ur_lon = w_ur_lon

    lbl02b = QLabel('latitude:')
    lbl02b.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl02b, irow, 2)

    w_ur_lat = QLineEdit()
    w_ur_lat.setFixedWidth(WDGT_SIZE_80)
    grid.addWidget(w_ur_lat, irow, 3)
    form.w_ur_lat = w_ur_lat

    # LL lon/lat
    # ==========
    irow += 1
    lbl01a = QLabel('Lower left longitude:')
    lbl01a.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl01a, irow, 0)

    w_ll_lon = QLineEdit()
    w_ll_lon.setFixedWidth(WDGT_SIZE_80)
    grid.addWidget(w_ll_lon, irow, 1)
    form.w_ll_lon = w_ll_lon

    lbl01b = QLabel('latitude:')
    lbl01b.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl01b, irow, 2)

    w_ll_lat = QLineEdit()
    grid.addWidget(w_ll_lat, irow, 3)
    w_ll_lat.setFixedWidth(WDGT_SIZE_80)
    form.w_ll_lat = w_ll_lat

    # bbox
    # ====
    irow += 1
    lbl03a = QLabel('AOI bounding box:')
    lbl03a.setAlignment(Qt.AlignRight)
    grid.addWidget(lbl03a, irow, 0)

    form.w_bbox = QLabel()
    grid.addWidget(form.w_bbox, irow, 1, 1, 5)

    return irow

def glblecss_limit_sims(form, grid, irow):
    """

    """
    irow += 1
    icol = 0
    lbl01 = QLabel('Limit simulations - start band, end band, max sims:')
    lbl01.setAlignment(Qt.AlignRight)
    helpText = 'limit number of latitude bands and/or number of simulation grids cells'
    lbl01.setToolTip(helpText)
    grid.addWidget(lbl01, irow, icol, 1, 2)

    icol += 2
    w_strt_band = QLineEdit()
    w_strt_band.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(w_strt_band, irow, icol)
    form.w_strt_band = w_strt_band

    icol += 1
    w_end_band = QLineEdit()
    w_end_band.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(w_end_band, irow, icol)
    form.w_end_band = w_end_band

    icol += 1
    w_max_sims = QLineEdit()
    w_max_sims.setFixedWidth(WDGT_SIZE_60)
    grid.addWidget(w_max_sims, irow, icol)
    form.w_max_sims = w_max_sims

    icol += 1
    w_view_run = QPushButton('View run')
    helpText = 'View details from last run'
    w_view_run.setToolTip(helpText)
    grid.addWidget(w_view_run, irow, icol)
    w_view_run.clicked.connect(form.viewRunReport)

    return irow

def grid_resolutions(form, grid, irow):
    """

    """
    irow += 1
    lbl16 = QLabel('Grid resolution:')
    lbl16.setAlignment(Qt.AlignRight)
    helpText = 'The size of each grid cell is described in arc minutes and arc seconds. The smallest cell resolution \n' \
        + 'corresponds to that of the HWSD database (30 arc seconds) and the largest to that used by the climate data ' \
        + '(30 arc minutes)'
    lbl16.setToolTip(helpText)
    grid.addWidget(lbl16, irow, 0)

    combo16 = QComboBox()
    for resol in sorted(RESOLUTIONS,reverse = True):
        combo16.addItem(str(RESOLUTIONS[resol]))
    combo16.setToolTip(helpText)
    combo16.setFixedWidth(WDGT_SIZE_60)
    form.combo16 = combo16
    combo16.currentIndexChanged[str].connect(form.resolutionChanged)
    form.combo16 = combo16
    grid.addWidget(combo16, irow, 1)

    form.lbl16a = QLabel('')
    form.lbl16a.setToolTip(helpText)
    grid.addWidget(form.lbl16a, irow, 2, 1, 3)

    return irow

def calculate_grid_cell(form, granularity = GRANULARITY):
    """

    """

    # use current lower left latitude for reference
    # =============================================
    latitude = 52.0     # default
    if hasattr(form, 'w_ll_lat'):
        latText = form.w_ll_lat.text()
        try:
            latitude = float(latText)
        except ValueError:
            print(latText)

    resol = form.combo16.currentText()
    try:
        granul = REVERSE_RESOLS[resol]
    except KeyError as err:
        print(str(err))
        return

    resol_deg = 1.0/float(granul)  # units of decimal degrees
    bbox = list([0.0, latitude, resol_deg, latitude + resol_deg])
    area = calculate_area(bbox)
    form.lbl16a.setText(format_bbox(bbox,area,2))

    form.req_resol_upscale = int(granularity/granul)    # number of base granuals making up one side of a cell
    form.req_resol_granul = granul                      # number of cells per degree
    form.req_resol_deg = resol_deg                      # size of each trapizoid cell in decimal degrees

    return resol_deg

def print_resource_locations(setup_file, config_dir, hwsd_dir, wthr_dir, lta_nc_fname, sims_dir, log_dir,
                             ecss_fns_dir=None):
    """
    report settings
    """
    print('\nResource locations:')
    print('\tsetup file:          ' + setup_file)
    print('\tconfiguration files: ' + config_dir)
    print('\tHWSD database:       ' + hwsd_dir)
    print('\tweather data sets:   ' + wthr_dir)

    if lta_nc_fname is not None:
        print('\tLTA weather file:    ' + lta_nc_fname)

    print('\tsimulations:         ' + sims_dir)
    print('\tlog_dir:             ' + log_dir)

    if ecss_fns_dir is not None:
        print('\tEcosse files:   ' + ecss_fns_dir)
    print('')

    return