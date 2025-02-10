#!/usr/bin/env python3

import argparse
import os
import os.path
import glob
import shutil
import time
import gzip
import math
import importlib.util
import pandas as pd

parser = argparse.ArgumentParser(
                    prog = 'merge_m2.py',
                    description = 'Takes m2 files from bam_to_m2.py aligned to same sequences and merges in order.',
                    epilog = 'Input/output are .*.txt.gz with comma-separated counts' )

parser.add_argument('-i','--infiles',nargs='*', required=True, help='names of *.txt.gz counts files to merge')
parser.add_argument('-o','--outfile',default='',help='name of output *.txt.gz file')
parser.add_argument('-v','--verbose',action="store_true",help='output as each sequence is processed')

args=parser.parse_args()

counts = []
df = None
df_init = False
numfiles = 0

infiles = []
for infile in args.infiles: infiles.extend(sorted(glob.glob(infile)))

time_readin = 0
time_output = 0
time_start = time.time()
time_startfile = time.time()
for (i,infile) in enumerate(infiles):

    # need to figure out separator
    sepchar = ','
    if infile.find('.gz')>-1: file_open=gzip.open
    else: file_open=open
    fid = file_open(infile,'rt')
    if fid.readline().find(sepchar) == -1: sepchar = ' '
    fid.close()

    # OK read it in!
    try:
        # uint32 needed to prevent memory overflow. Note maximum is 4Gb, so this will soon be issue, and need to shift to streaming.
        if args.verbose: print('Reading: ',infile)
        df_infile = pd.read_table(infile,sep=sepchar,header=None,dtype='uint32')
        df_infile_read_correctly = True
    except pd.errors.EmptyDataError:
        print('Note: %s was empty. Skipping.' % infile )
        df_infile_read_correctly = False
        continue # will skip the rest of the block and move to next file
    if not df_init and df_infile_read_correctly:
        df = df_infile
        df_init = True
    elif df_infile_read_correctly:
        df = df.add(df_infile,fill_value = 0)
    numfiles += 1

time_after_infile = time.time()

time_readin += time_after_infile - time_startfile

if df is not None:
    outfile = args.outfile
    if args.verbose: print('Writing: ',outfile)
    compression=None
    if len(outfile)>3 and outfile[-3:]=='.gz': compression='gzip'
    df.to_csv(outfile,sep=',',header=None,index=None,compression=compression)

    nseq = len(df)
    tot_counts = df.max(axis=1).sum()
    print( 'Compiled %8d total counts for %6d sequences from %6d of %6d files into: %s' % (tot_counts,nseq,numfiles,len(infiles),outfile) )

time_output += (time.time()-time_after_infile)

time_end = time.time()


print( '\nTimings:')
print( 'Read in data: ' + time.strftime("%H:%M:%S",time.gmtime(time_readin) ) )
print( 'Output  data: ' + time.strftime("%H:%M:%S",time.gmtime(time_output) ) )

print( '\nTotal time: ' + time.strftime("%H:%M:%S",time.gmtime(time_end-time_start) ) )

