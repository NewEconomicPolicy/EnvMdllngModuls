#-------------------------------------------------------------------------------
# Name:        jon_fns.py
# Purpose:     script to create objects describing NC data sets
# Author:      Mike Martin
# Created:     31/05/2020
# Licence:     <your licence>
#-------------------------------------------------------------------------------
#!/usr/bin/env python

__prog__ = 'jon_fns.py'
__version__ = '0.0.0'

# Version history
# ---------------
#
from shutil import copytree
from os import listdir, makedirs, mkdir
from os.path import exists, split, join, isdir
from time import time
from datetime import timedelta
from win32api import GetVolumeInformation

from PyQt5.QtWidgets import QApplication

ERROR_STR = '*** Error *** '
WARNING_STR = '*** Warning *** '
N_DECIM = 4
sleepTime = 5

def copy_jon_lta_data(form, use_drive, out_drive):
    """
    C
    """
    lta_inp_dir = join(use_drive, 'ECOSSE_LTA')
    lta_out_dir = join(out_drive, 'PortableSSD', 'ECOSSE_LTA')

    if not isdir(lta_out_dir):
        makedirs(lta_out_dir)

    icntr = 0
    for rcp in listdir(lta_inp_dir):
        dirnm_inp = join(lta_inp_dir, rcp)
        dirnm_out = join(lta_out_dir, rcp)
        print('Copying lta dir: ' + dirnm_inp + ' to ' + dirnm_out)
        t1 = time()
        if exists(dirnm_out):
            print('Output directory: ' + dirnm_out + ' already exists')
            continue
        copytree(dirnm_inp, dirnm_out)
        t2 = time()
        scnds_elapsed = round(t2 - t1)
        icntr += 1
        print('Finished copying LTA ' + rcp + '\tTime taken: ' + str(timedelta(seconds=scnds_elapsed)))
        if icntr > 0:
            break

def copy_jon_wthr_data(form, use_drive, out_drive):
    """
     assumption if that SSD data is consistent
    """
    wthr_inp_dir = join(use_drive, 'ECOSSE_RCP')
    wthr_out_dir = join(out_drive, 'PortableSSD', 'ECOSSE_RCP')
    max_recs = form.w_max_recs.text()
    max_recs = int(max_recs)

    last_time = time()
    strt_time = time()
    ncopied_all = 0
    for rcp in listdir(wthr_inp_dir):
        dirnm_inp = join(wthr_inp_dir, rcp)
        dirs_to_copy = listdir(dirnm_inp)
        ndirs2cpy = len(dirs_to_copy)

        dirnm_out = join(wthr_out_dir, rcp)
        print('Copying wthr dir: ' + dirnm_inp + ' to ' + dirnm_out)
        ncopied = 0
        for coord in dirs_to_copy:
            coord_dir_inp = join(dirnm_inp, coord)
            coord_dir_out = join(dirnm_out, coord)
            if exists(coord_dir_out):
                continue
            else:
                last_time = _update_progress(last_time, form.w_prgrss, rcp, ncopied, ndirs2cpy)
                copytree(coord_dir_inp, coord_dir_out)

                ncopied += 1
                ncopied_all += 1
                if ncopied >= max_recs:
                    break

    scnds_elapsed = round(time() - strt_time)
    mess = '\nFinished after N copies: {}\ttime taken: '.format(ncopied_all)
    mess += str(timedelta(seconds=scnds_elapsed))
    form.w_prgrss.setText(mess)

    return

def _update_progress(last_time, w_prgrss, rcp, ncopied, ndirs2cpy):
    """

    """
    new_time = time()
    if new_time - last_time > sleepTime:
        prcnt_cells = round(100* (ncopied/ndirs2cpy), 2)
        mess = 'RCP: {}\tCopied: {}  cells\t% processed: {}'.format(rcp, ncopied, prcnt_cells)
        w_prgrss.setText(mess)
        QApplication.processEvents()
        last_time = new_time

    return last_time

def create_bash_script(form, san_disk_drv, out_drv = 'F:\\'):
    """
    write a linux script to copy data from SanDisk
    """
    wthr_inp_dir = join(san_disk_drv, 'ECOSSE_RCP')
    from_drv = san_disk_drv.lower()[0]
    to_drv = out_drv.lower()[0]

    print('\n')

    out_recs = []

    for dn in ['ECOSSE_LTA', 'ECOSSE_RCP', 'temp']:
        out_dir = join(out_drv, 'PortableSSD', dn)
        if not isdir(out_dir):
            makedirs(out_dir)

    for rcp in listdir(wthr_inp_dir):
        if rcp.find('rcp') == 0:
            cmd_str = 'cp -pr /mnt/' + from_drv + '/ECOSSE_LTA/' + rcp + ' /mnt/' + to_drv + '/PortableSSD/ECOSSE_LTA'
            out_recs.append(cmd_str + '\n')

    out_recs.append('\n')

    for rcp in listdir(wthr_inp_dir):
        if rcp.find('rcp') == 0:
            cmd_str = 'cp -pr /mnt/' + from_drv + '/ECOSSE_RCP/' + rcp + ' /mnt/' + to_drv + '/PortableSSD/ECOSSE_RCP'
            out_recs.append(cmd_str + '\n')

    fn = join(out_dir, 'bash_script.sh')
    with open(fn, 'w') as fbash:
        fbash.writelines(out_recs)
    print('Wrote:  ' + fn)

    return

def identify_ssd(prtbl_ssd = 'PortableSSD'):
    '''
    check SSD is accessible
    '''
    ssd_found = False
    for drive in range(ord('D'), ord('Z')):
        drv_chr = chr(drive)
        if exists(drv_chr + ':'):
            drive_win = drv_chr + ':\\'
            vol_info = GetVolumeInformation(drive_win)
            if vol_info[0] == prtbl_ssd:
                use_drive = drive_win
                ssd_found = True
                break

    return ssd_found, use_drive, prtbl_ssd
