"""
#-------------------------------------------------------------------------------
# Name:
# Purpose:     create masks for pasture (33) and cropland (22)
# Author:      Mike Martin
# Created:     11/12/2015
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#
"""

__prog__ = 'hilda_fns.py'
__version__ = '0.0.1'
__author__ = 's03mm5'

from time import time
from os.path import split, join, isfile
from os import remove
from netCDF4 import Dataset, stringtochar
from numpy import arange, float64, array
from math import floor, ceil

from mask_fns_and_class import LanduseSet, make_mask_nc, fetch_prevalence
from country_fns import make_country_nc, update_progress

N_DECIM = 4
MAX_CELLS = 9999999999

SOIL_FN = 'soildata_4000.txt'
SOIL_FN = 'soildata.txt'

def add_codes_countries_nc(form):
    '''
    called from GUI - read Astley file:
    '''
    glec_dir, dummy = split(form.hwsd_dir)
    soil_dir = join(glec_dir, 'Soil')
    states_fn = join(soil_dir, 'all_Countries.nc')
    if not isfile(states_fn):
        print('Countries NC file ' + states_fn + ' must exist')
        return

    acountries, cnames = soildata_codes_vault()
    ncntries = len(acountries)

    nc_dset = Dataset(states_fn, mode='a', format='NETCDF4')

    # create country integer code variable
    # ====================================
    nc_dset.createDimension('cntry_code', ncntries)

    codes = nc_dset.createVariable('cntry_code', 'i2', ('cntry_code',))
    codes.description = 'country codes as used in miscanfor.f90'
    codes.long_name = 'Country codes'
    codes[:] = acountries

    # ascertain max length of country names
    # =====================================
    maxlen = 0
    for cname in cnames:
        maxlen = max(maxlen, len(cname))
    spacer = maxlen*' '
    spcr_fmt = 'S' + str(maxlen)

    # create country names variable
    # =============================
    nc_dset.createDimension('cntry_name', ncntries)

    names = nc_dset.createVariable('cntry_name', spcr_fmt, ('cntry_name',))
    names.description = 'country names as used in miscanfor.f90'
    names.long_name = 'Country names'

    # use stringtochar function to convert a string array to a character array
    # ========================================================================
    char_cnames = [stringtochar(array([cname], spcr_fmt)) for cname in cnames]
    # names[:] = char_cnames
    names[:] = cnames

    nc_dset.close()

    return

def write_countries_nc(form):
    '''
    called from GUI - read Astley file:
    '''
    glec_dir, dummy = split(form.hwsd_dir)
    soil_dir = join(glec_dir, 'Soil')
    csv_fname = join(soil_dir, SOIL_FN)
    if not isfile(csv_fname):
        mess = 'CSV file of states ' + csv_fname + ' must exist... \n'
        print(mess)
        form.lgr.info(mess)
        return

    fobj = open(csv_fname, 'r')
    resol_rec = fobj.readline()
    fobj.close()
    dummy, resol_str = resol_rec.split()
    resol = float(resol_str)

    # countries NC file
    # =================
    states_fn = join(soil_dir, 'all_Countries.nc')
    if isfile(states_fn):
        try:
            remove(states_fn)
        except PermissionError as err:
            print(str(err))
            return None

        print('Deleted ' + states_fn)

    make_country_nc(csv_fname, states_fn, resol)

    return

def _test_area(resol_d2, resol, max_cells):
    '''
    generate european subset
    '''
    lon_ll = 15.9
    lat_ll = 44.9
    lon_ur = 27.32
    lat_ur = 52.7

    alons = arange(floor(lon_ll) + resol_d2, ceil(lon_ur) - resol_d2, resol, dtype=float64)
    lons = [round(float(lon), N_DECIM) for lon in alons]

    alats = arange(floor(lat_ll) + resol_d2, ceil(lat_ur) - resol_d2, resol, dtype=float64)
    lats = [round(float(lat), N_DECIM) for lat in alats]

    mess = 'Will generate maximum of {} cells:\n\t'.format(max_cells)
    mess += 'Lower lat/long: {} {}\t\tUpper lat/long: {} {}'.format(lat_ll, lon_ll, lat_ur, lon_ur)
    print(mess)

    return lats, lons

def generate_masks(form):
    '''
    called from GUI
    '''
    wthr_rsce = form.combo10w.currentText()
    wthr_rsce_key  = wthr_rsce + '_Mnth'
    wthr_set = form.weather_sets[wthr_rsce_key]
    
    # make sure bounding box is correctly set
    lon_ll = wthr_set['lon_ll']
    lat_ll = wthr_set['lat_ll']
    lon_ur = wthr_set['lon_ur']
    lat_ur = wthr_set['lat_ur']
    bbox_aoi =  list([lon_ll, lat_ll, lon_ur, lat_ur])

    resol = form.req_resol_deg
    resol_d2 = resol/2

    # open land use NC file
    # =====================
    glec_dir, dummy  = split(form.hwsd_dir)

    hilda_dir = join(glec_dir, 'Hilda_land_use', 'hildap_vGLOB-1.0-f')

    lu_fname = join(hilda_dir, 'hildaplus_vGLOB-1.0-f_states.nc')
    if not isfile(lu_fname):
        mess = 'HILDA land use file ' + lu_fname + ' must exist... \n'
        print(mess); form.lgr.info(mess)
        return

    lu_defn = LanduseSet(lu_fname)
    lu_defn.nc_dset = Dataset(lu_fname, mode='r')

    mask_fn = make_mask_nc(hilda_dir, resol, lon_ll, lat_ll, lon_ur, lat_ur)
    if mask_fn is None:
        return
    mask_defn = LanduseSet(mask_fn)
    mask_defn.nc_dset = Dataset(mask_fn, mode='a', format='NETCDF4')

    # for each lat/lon derive prevalence of land use 22 and 33
    # ========================================================
    form.lgr.info('\n\n')
    lats = mask_defn.lats
    lons = mask_defn.lons

    # lats, lons = _test_area(resol_d2, resol, MAX_CELLS)

    nvoids = 0
    nerrors = 0
    no_squares = 0
    last_time = time()
    start_time = time()
    completed = 0
    ncells = len(lats)*len(lons)

    for lat in lats:
        for lon in lons:
            fetch_prevalence(form.lgr, lu_defn, mask_defn, lat, lon, resol_d2, nvoids, no_squares, nerrors)
            completed += 1
            if completed >= MAX_CELLS:
                break
            last_time = update_progress(last_time, start_time, completed, ncells, nvoids, no_squares)

        if completed >= MAX_CELLS:
            print('\nCompleted after creating {} cells'.format(completed))
            break

    # close NC files
    # ==============
    mask_defn.nc_dset.close()
    lu_defn.nc_dset.close()

    return

def soildata_codes_vault():
    '''
    from miscanfor.f90
    '''
    acountries = [58, 32, 11, 218, 165, 183, 232, 107, 54, 173, 200, 22, 101, 75, 72, 133, 70, 39, 18, 112, 17, 64, 1,
                  31, 33, 143, 66, 38, 13, 42, 160, 151, 241, 242, 26, 25, 87, 216, 117, 188, 153, 152, 139, 12, 123,
                  115, 71, 23, 211, 53, 68, 78, 222, 138, 167, 113, 150, 240, 111, 239, 163, 210, 213, 88, 230, 5, 149,
                  7, 52, 110, 14, 99, 77, 102, 170, 0, 145, 144, 142, 76, 114, 27, 80, 100, 221, 129, 51, 79, 225, 47,
                  234, 84, 83, 192, 48, 55, 157, 82, 136, 132, 61, 65, 36, 49, 28, 10, 21, 219, 35, 19, 215, 34, 105,
                  46, 98, 9, 122, 171, 120, 187, 103, 41, 217, 63, 8, 154, 134, 164, 62, 174, 176, 229, 202, 116, 16,
                  147, 109, 130, 59, 166, 140, 137, 2, 141, 89, 108, 104, 168, 74, 124, 97, 43, 69, 207, 172, 20, 148,
                  127, 119, 146, 73, 29, 24, 94, 156, 106, 81, 85, 126, 4, 45, 86, 220, 50, 92, 60, 155, 93, 15, 91,
                  121, 90, 56, 205, 44, 125, 131, 226, 228, 3, 57, 161, 175, 135, 236, 6, 128, 30, 37, 40]

    cnames = ['Afghanistan', 'Albania', 'Algeria', 'Andorra', 'Angola', 'Anguilla', 'Antigua and Barbuda',
              'Argentina'] + \
             ['Armenia', 'Aruba', 'Australia', 'Austria', 'Azerbijan', 'Bahamas', 'Bahrain', 'Bangladesh', 'Barbados',
              'Belarus', 'Belgium', 'Belize'] + \
             ['Benin', 'Bhutan', 'Bolivia', 'Bosnia and Herezegovina', 'Botswana', 'Brazil', 'Brunei', 'Bulgaria',
              'Burkina Faso', 'Burundi', 'Cambodia', 'Cameroon'] + \
             ['Canada', 'Canada', 'Central African Repbulic', 'Chad', 'Chile', 'China', 'Colombia', 'Comoros',
              'Congo (Democratic Republic)', 'Congo (Republic)'] + \
             ['Costa Rica', 'Cote D\'Ivoir', 'Croatia', 'Cuba', 'Cyprus', 'Czech Republic', 'Denmark', 'Djibouti',
              'Dominica', 'Dominican Republic', 'East Timor', 'Ecuador'] + \
             ['Egypt', 'El Salvador', 'Equatorial Guinea', 'Eritrea', 'Estonia', 'Ethiopia', 'Falkland Islands',
              'Faroe Islands', 'Fiji', 'Finland', 'France', 'French Guiana'] + \
             ['Gabon', 'Gambia', 'Georgia', 'Germany', 'Ghana', 'Greece', 'Greenland', 'Grenada', 'Guadeloupe',
              'Guatemala', 'Guinea', 'Guinea-Bissau', 'Guyana', 'Haiti', 'Honduras'] + \
             ['Hungary', 'Iceland', 'India', 'Indonesia', 'Iran', 'Iraq', 'Ireland', 'Isle of Man', 'Israel', 'Italy',
              'Jamaica', 'Japan', 'Jersey', 'Jordan', 'Kazakhstan', 'Kenya'] + \
             ['Korea (North)', 'Korea (South)', 'Kuwait', 'Kyrgyzstan', 'Laos', 'Latvia', 'Lebanon', 'Lesotho',
              'Liberia', 'Libya', 'Lichtenstein', 'Lithuania', 'Luxembourg'] + \
             ['Macau', 'Macedonia', 'Madagascar', 'Malawi', 'Malaysia', 'Mali', 'Malta', 'Martinique', 'Mauritania',
              'Mayotte', 'Mexico', 'Moldova', 'Monaco', 'Mongolia', 'Morocco'] + \
             ['Mozambique', 'Myanmar', 'Namibia', 'Nepal', 'Netherlands', 'Netherlands Antilles', 'New Caledonia',
              'New Zealand', 'Nicaragua', 'Niger', 'Nigeria', 'Norway', 'Oman'] + \
             ['Pakistan', 'Palestine', 'Panama', 'Papua New Guinea', 'Paraguay', 'Peru', 'Philippines', 'Poland',
              'Portugal', 'Puerto Rico', 'Qatar', 'Romania', 'Russia', 'Rwanada'] + \
             ['Saint Lucia', 'Saint Pierre and Miquelon', 'Saint Vincent and the Grenadines', 'San Marino',
              'Sao Tome and Principe', 'Saudi Arabia', 'Senegal'] + \
             ['Sierra Leone', 'Singapore', 'Slovakia', 'Slovenia', 'Soloman Islands', 'Somalia', 'South Africa',
              'Spain', 'Sri Lanka', 'Sudan', 'Suriname', 'Swaziland'] + \
             ['Sweden', 'Switzerland', 'Syria', 'Taiwan', 'Tajikistan', 'Tanzania', 'Thailand', 'Togo',
              'Trinidad and Tobago', 'Tunisia', 'Turkey', 'Turkmenistan'] + \
             ['Turks and Caicos Islands', 'Uganda', 'Ukraine', 'United Arab Emirates', 'United Kingdom',
              'United States', 'Uruguay', 'Uzbekistan', 'Vanuatu'] + \
             ['Venezuela', 'Vietnam', 'Virgin Islands (U.S.)', 'Western Sahara', 'Yemen', 'Yugoslavia', 'Zambia',
              'Zimbabwe']

    return acountries, cnames