from shutil import copytree
from os import remove, rmdir, listdir
from time import time
from os.path import join, isdir, isfile, splitdrive, splitext
t1 = time()
strt_time = t1
iskip, icopy = 2*[0]
for shrt_dir in drns_in:
    dir_in = join(rcp_dir_inp, shrt_dir)
    dir_out = join(rcp_dir_out, shrt_dir)
    '''
    if isdir(dir_out):
        remove(dir_out) 
    '''
    try:
        copytree(dir_in, dir_out)
        icopy += 1
    except FileExistsError as err:
        iskip += 1
    
    t2 = time()    
    lapsed_time = t2 - t1
    if lapsed_time > 8:
        print('copied {}\tskipped: {}\tcoords'.format(icopy, iskip))
        t1 = t2
    
print('copied {}\tskipped: {}\tcoords in {} secs'.format(icopy, iskip, int(time()-strt_time)))