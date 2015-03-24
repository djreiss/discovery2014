import os
import glob
import numpy as np
import pandas as pd

import utils as ut

def parse_args():
    from optparse import OptionParser
    parser = OptionParser()
    ## None = process all samples in remote_dir
    parser.add_option('-s', '--sample', dest='sample',
                        help='sample to process', metavar='SAMPLE', default=None) #'Sample_US-1505885')
    parser.add_option('-r', '--remote_dir', dest='remote_dir', default='./ALL_SAMPLES/', 
                        help='remote directory [default: %default]', metavar='REMOTE_DIR')
    parser.add_option('--pylab') ## ignore this if run from ipython
    
    options, args = parser.parse_args()
    return(options, args)

try: ## fails in ipython if called via 'import params'
    options,args = parse_args()
    print options.sample
    print options.remote_dir
except:
    None

sample = None
try:
    sample = options.sample
except:
    sample = None

## directory that remote sshfs directory with samples is mounted as
remote_dir = './ALL_SAMPLES/'
try:
    remote_dir = options.remote_dir + '/'
except:
    remote_dir = './ALL_SAMPLES/'

genome_lookup = {
    'Desulfovibrio_vulgaris_Hildenborough_uid57645':['NC_002937','NC_005863'],
    'Desulfovibrio_alaskensis_G20_uid57941':['NC_007519'],
    'Methanococcus_maripaludis_S2_uid58035':['NC_005791'],
    'Methanospirillum_hungatei_JF_1_uid58181':['NC_007796'] }

sample_tab = pd.read_excel( 'ListAllSamplesRNAseq_Discovery_2014.xlsx', 1 )
orgs = sample_tab.organisms.value_counts()
org_lookup = {
    'D vulgaris Hildenborugh':['Desulfovibrio_vulgaris_Hildenborough_uid57645'],
    'D alaskensis G20':['Desulfovibrio_alaskensis_G20_uid57941'],
    'D vulgaris Hildenborugh and M maripaludis': \
            ['Desulfovibrio_vulgaris_Hildenborough_uid57645','Methanococcus_maripaludis_S2_uid58035'],
    'D alskensis G20 and M maripaludis': \
            ['Desulfovibrio_alaskensis_G20_uid57941','Methanococcus_maripaludis_S2_uid58035'],
    'D alskensis G20 and Methanospirillium hangeitii': \
            ['Desulfovibrio_alaskensis_G20_uid57941','Methanospirillum_hungatei_JF_1_uid58181'],
    'D vulgaris Hildenborugh and Methanospirillium hangeitii': \
            ['Desulfovibrio_vulgaris_Hildenborough_uid57645','Methanospirillum_hungatei_JF_1_uid58181'] }

sample_files = np.sort( np.array( glob.glob( remote_dir + 'dvh-output/Sample*' ) ) )
print sample_files
if sample is not None:
    sample_files = sample_files[ np.where( [os.path.basename(f)==sample for f in sample_files] ) ]
print sample_files

## also read in the Sample*/Log.final.out file for statistics

all_counts = {}
logs = {}
for org in genome_lookup.keys():
    all_counts[org] = {}
    logs[org] = {}

    for SAMPLE in sample_files:
        #set SAMPLE="Sample_US-1505885"

        SAMP = os.path.basename(SAMPLE)
        #all_counts[org][SAMP] = {}
        #logs[org][SAMP] = {}

        counts_files = glob.glob( SAMPLE + '/' + org + '/*.counts' )
        if len(counts_files) <= 0:
            continue
    
        xx = {}
        log = None
        for counts_file in counts_files:
            print counts_file
            x = pd.read_table(counts_file, header=None) ## index_col=0, 
            print x.sum(0)[1]
            xx[os.path.basename(counts_file).replace('.counts','')] = x ##pd.concat( [xx, x], axis=0 )

            log_file = os.path.dirname(counts_file) + '/Log.final.out'
            log = pd.read_table( log_file, index_col=0, header=None ) ##ut.readLines( log_file )

        xx = pd.concat( xx )
        xx = xx.set_index( xx[0] ); xx = xx.drop([0], 1)
        all_counts[org][SAMP] = xx
        logs[org][SAMP] = log

all_counts = {i:pd.concat(j,1) for i,j in all_counts.items()}
logs = {i:pd.concat(j,1) for i,j in logs.items()}

writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_counts.xlsx')
for k in all_counts.keys():
    all_counts[k].to_excel( writer, k[0:30] )
writer.save()

writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_logs.xlsx')
for k in logs.keys():
    logs[k].to_excel( writer, k[0:30] )
writer.save()

