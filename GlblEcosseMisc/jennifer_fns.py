#-------------------------------------------------------------------------------
# Name:        filter_hwsd_fns
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'jennifer_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from pandas import read_csv, DataFrame
from numpy import isnan, nan
from os import remove
from os.path import exists, split
from time import time

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
N_DECIM = 4

def filter_openfoodfacts_csv(form):
    """

    """
    max_recs = int(form.w_max_recs.text())

    dset_inp_fn = 'E:\\Jennifer\\en_food_1000.csv'
    dset_inp_fn = 'E:\\Jennifer\\en_openfoodfacts_org_products.csv'
    fn_vld_indcs =  'E:\\Jennifer\\en_openfoodfacts_valid_recs.csv'
    dset_out_fn = 'E:\\Jennifer\\en_openfoodfacts_summary.csv'

    col_list = list(range(0, 200))
    rqrd_cols = (col_list[:2] + [col_list[8]] + [col_list[11]] + col_list[13:16] + [col_list[17]] +
                 col_list[19:21] + col_list[22:24] + [col_list[25]] + col_list[27:29])

    rqrd_cols += (col_list[34:36] + col_list[38:42] + [col_list[50]] + col_list[52:57] + col_list[59:62] +
                  col_list[65:70] + col_list[77:79] + col_list[85:91])

    rqrd_cols += (col_list[126:129] + [col_list[138]] + [col_list[141]] + col_list[145:148] + col_list[183:187] +
                  col_list[190:195])

    '''
    print('Counting lines in dataset: ' + dset_inp_fn)
    with open(dset_inp_fn) as fobj:
        for ic, dummy in enumerate(fobj.readline()):
            pass
    '''

    # df_valid = DataFrame(dset_inp_fn, columns = rqrd_cols)
    # df_out = read_csv(dset_inp_fn, sep='\t', usecols = rqrd_cols, nrows=0)

    strt_time = time()
    print('Loading {:,} lines from dataset {}, this can take some time...'.format(max_recs, dset_inp_fn))
    df = read_csv(dset_inp_fn, nrows = max_recs, sep='\t', usecols=rqrd_cols)
    t_load = round((time() - strt_time) / 60, 2)
    print('Time to load ' + split(dset_inp_fn)[1] + ' {} mins'.format(t_load))

    # stage 1 - create list of valid records
    # ======================================

    del_indices = []
    valid_indices = []
    indx = 0

    # check rows 56, 57 and 68
    # ========================
    for nutri, nova, eco in zip(df['nutriscore_grade'].values, df['nova_group'].values, df['ecoscore_grade'].values):

        if type(nutri) is str:
            nutri = 999  # indicates that nutri will be interpreted as string

        if isnan(nutri) or isnan(nova) or eco == 'unknown':
            # print('No data for {}'.format(indx))
            del_indices.append(indx)
        else:
            valid_indices.append(indx)

        indx += 1
        if 100000 * int(indx / 100000) == indx:
            print('Found {} records with no data\t{} with data\tindx: {}'.format(len(del_indices),
                                                                                 len(valid_indices), indx))

        if indx > max_recs:
            break

    nvalids = len(valid_indices)
    print('*** Finished processing file {}\n\tfound {:,} records with no data\t{:,} with data\n'.format(dset_inp_fn,
                                                                     len(del_indices), nvalids))

    df_out = df.iloc[valid_indices]

    print('*** Finished stage 2 after processing {:,} records\n'.format(nvalids))

    print('Writing ' + dset_out_fn)
    if exists(dset_out_fn):
        remove(dset_out_fn)

    df_out.to_csv(dset_out_fn, sep = '\t', index=True, header = True)

    print('*** Filtering finished ***')

    return

def edit_mngmnt(form):
    """
    C
    """
    fn = 'E:\\temp\\superg\\management.txt'
    with open(fn) as fobj:
        lines = fobj.readlines()

    new_lines = []
    for line in lines:
        if line.find('Amount of manure') > 0:
            new_lines.append('5.0  ' + line[5:])
        else:
            new_lines.append(line)

    fn_out = 'E:\\temp\\superg\\management_out.txt'
    with open(fn_out, 'w') as fobj:
        fobj.writelines(new_lines)
