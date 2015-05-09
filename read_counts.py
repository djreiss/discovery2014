import os
import glob
import numpy as np
import pandas as pd

import utils as ut

if not 'DO_SAVE' in globals():
    DO_SAVE = False

print DO_SAVE

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

if DO_SAVE:
    writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_counts.xlsx')
    for k in all_counts.keys():
        all_counts[k].to_excel( writer, k[0:30] )
    writer.save()

    writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_logs.xlsx')
    for k in logs.keys():
        logs[k].to_excel( writer, k[0:30] )
    writer.save()



## Now normalize the counts - first by millions of (total) transcripts, then by length_of_gene/mean_length

#all_freqs = { i: j.copy() / j.sum(0) * 1.0e6 for i,j in all_counts.items() }

## Don't use read counts that were not assigned to anything (for totals).
## Don't use read counts that were assigned to rRNAs (for totals).

all_freqs = {}

import gff3  ## Get length of gene from gff files
for org in all_counts.keys():
    print 'NORMALIZING: ', org
    freq_mat = all_counts[ org ].copy()
    gene_lengths = pd.Series()
    is_rna = pd.Series()   ## note this is only rRNAs -- not tRNAs (don't have as large copy numbers)
    for i in genome_lookup[ org ]:
        gff = 'GENOMES/' + i + '.gff'
        ##tab = pd.read_table( gff, index_col=0, sep='\t', comment='#', header=None )
        ##recs = [ gff3.parseGFFAttributes( tab.ix[i][8] ) for i in range(1,tab.shape[0]) ]
        recs = { rec.attributes['ID']: rec for rec in gff3.parseGFF3( gff ) }
        for geneid,rec in recs.items():
            isRNA = rec.type == 'rRNA' ## or rec.type == 'transcript' or rec.type == 'tRNA' 
            while not rec.attributes.has_key( 'locus_tag' ) and rec.attributes.has_key('Parent'):
                geneid,rec = rec.attributes['Parent'], recs[rec.attributes['Parent']]
            if not rec.attributes.has_key( 'locus_tag' ):
                continue
            if isRNA or not rec.attributes[ 'locus_tag' ] in is_rna.index:
                is_rna[ rec.attributes[ 'locus_tag' ] ] = isRNA
            gene_lengths[ rec.attributes[ 'locus_tag' ] ] = np.float( np.abs( rec.end - rec.start ) )
                
    gene_lengths = gene_lengths[ freq_mat.index ].fillna( 1 ) ## also sets e.g. '__no_feature' to 1
    sum_reads = freq_mat.sum().sum()
    sum_aligned = freq_mat.ix[is_rna.index].sum().sum()
    is_rna = is_rna[ freq_mat.index ] ## also sets e.g. '__no_feature' to is_rna=True
    sum_rrna = freq_mat.ix[is_rna == True].sum().sum()
    is_rna = is_rna.fillna( True ) ## also sets e.g. '__no_feature' to is_rna=True
    is_valid = (is_rna == True) | (is_rna == False)
    print is_rna.sum(), 'rRNAs', sum_reads, sum_aligned, sum_rrna
    mean_len = gene_lengths[ ~is_rna ].mean()
    freq_mat = ( freq_mat / freq_mat[ ~is_rna ].sum(0) * 1.0e6 ).div( gene_lengths, axis='index' ) * mean_len
    all_freqs[org] = freq_mat

if DO_SAVE:
    writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_freqs.xlsx')
    for k in all_freqs.keys():
        np.round(all_freqs[k], 2).to_excel( writer, k[0:30] )
    writer.save()





## Now get the sample info

sample_info = pd.read_excel('Sample_Info_COMPLETE.xlsx') ##,skiprows=[0])
sample_info = sample_info.set_index( sample_info['Sample name'] )
sample_infos = { k:sample_info.ix[all_freqs[k].columns.droplevel(1).values] for k in all_freqs.keys() }

if DO_SAVE:
    writer = pd.ExcelWriter(remote_dir + 'dvh-output/all_sample_infos.xlsx')
    for k in sample_infos.keys():
        sample_infos[k].to_excel( writer, k[0:30] )
    writer.save()



raise 'STOPPING'


## Try plotting an operon - DVU000[2-5]
from matplotlib import pyplot as plt
pd.options.display.mpl_style = 'default'
from mpltools import style
from mpltools import layout
style.use('ggplot')

#x = all_counts['Desulfovibrio_vulgaris_Hildenborough_uid57645'].copy()

x = all_freqs['Desulfovibrio_vulgaris_Hildenborough_uid57645'].copy()
x[ x == 0 ] = np.NAN
x = np.log10( x )
x = (x.transpose() - x.mean(1)).transpose()
x = (x.transpose() / x.std(1)).transpose()
## another option: sklearn.preprocessing.scale(x.values, axis=0, with_mean=True, with_std=True, copy=True)
op = x.ix[ ['DVU0002','DVU0003','DVU0004','DVU0005','DVU0018'] ] 
#op = x.ix[ ['DVU0060','DVU0061'] ] 
op.transpose().plot(kind='line')

## samples where values are dropping way low:
samps = op.columns.droplevel(1).values[np.where((op < 1.5).sum()>0)]
sample_infos['Desulfovibrio_vulgaris_Hildenborough_uid57645'].ix[ samps ] ## all starvation...

from sklearn.cluster import KMeans
inertias = np.zeros(400)
for n_clust in range(50,400):
    km = KMeans(init='k-means++', n_clusters=n_clust, n_init=200, n_jobs=6)
    ## k-means cant handle nans so replace with zeros
    km.fit( x.fillna( 0 ).values )
    #op = x.ix[ km.labels_ == 0 ]
    #op.transpose().plot(kind='line')
    inertias[n_clust] = km.inertia_
    print n_clust, km.inertia_
pd.Series(inertias).plot()


