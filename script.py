import os
import glob
import tempfile
import numpy as np
import utils as ut

def parse_args():
    from optparse import OptionParser
    parser = OptionParser()
    ## None = process all samples in remote_dir
    parser.add_option('-s', '--sample', dest='sample',
                        help='sample to process', metavar='SAMPLE', default=None) #'Sample_US-1505885')
    parser.add_option('-r', '--remote_dir', dest='remote_dir', default='./ALL_SAMPLES/', 
                        help='remote directory [default: %default]', metavar='REMOTE_DIR')
    parser.add_option('-o', '--organism', dest='organism', default=None, 
                        help='organism specifications, comma separated [default: %default]', metavar='ORGANISMS')
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

## fname here is something like 'Desulfovibrio_vulgaris_Hildenborough_uid57645/NC_002937'
## -- prefix of fna file from ftp://ftp.ncbi.nih.gov/genomes/Bacteria/
#ftp://ftp.ncbi.nih.gov/genomes/Bacteria/Desulfovibrio_alaskensis_G20_uid57941/NC_007591

## TBD: download all genomes/plasmids for a given organism instead of having to specify each one

def dld_genome_data( fname, url_prefix='ftp://ftp.ncbi.nih.gov/genomes/Bacteria/' ):
    bname = os.path.basename( fname ) ## e.g. NC_002937

    try:
        os.makedirs( 'GENOMES' )
    except:
        pass

    if not os.path.isfile( 'GENOMES/' + bname + '.gff' ):
        try:
            os.remove('GENOMES/' + bname + '.gff'); os.remove('GENOMES/' + bname + '.fna'); 
        except:
            pass

        os.chdir( 'GENOMES' )
        cmd = "wget '" + url_prefix + fname + ".gff'" 
        tmp = ut.system( cmd )
        cmd = "wget '" + url_prefix + fname + ".fna'"
        tmp = ut.system( cmd )

        myfile = open( bname + '.fna', 'r' )
        line = myfile.readline()
        myfile.close()
        cmd = "../rpl '" + line.strip() + "' '>" + line.split('|')[3] + "' " + bname + '.fna'
        tmp = ut.system( cmd )

        ut.system( "awk '($3==\"gene\"){print}' " + bname + ".gff > " + bname + "_genes.gff" )
        os.chdir( '..' )

        #rpl 'ID=' 'gene_id=' NC_002937_genes.gff
        ## not necessary?:
        #ut.system( "cat NC_002937_genes.gff NC_005863_genes.gff >DvH_genes.gff" )

genome_lookup = {
    'Desulfovibrio_vulgaris_Hildenborough_uid57645':['NC_002937','NC_005863'],
    'Desulfovibrio_alaskensis_G20_uid57941':['NC_007519'],
    'Methanococcus_maripaludis_S2_uid58035':['NC_005791'],
    'Methanospirillum_hungatei_JF_1_uid58181':['NC_007796'] }

dld_genome_data( 'Desulfovibrio_vulgaris_Hildenborough_uid57645/NC_002937' )
dld_genome_data( 'Desulfovibrio_vulgaris_Hildenborough_uid57645/NC_005863' )
dld_genome_data( 'Desulfovibrio_alaskensis_G20_uid57941/NC_007519' )
dld_genome_data( 'Methanococcus_maripaludis_S2_uid58035/NC_005791' )
dld_genome_data( 'Methanospirillum_hungatei_JF_1_uid58181/NC_007796' )

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

gff_files = np.array( os.listdir( 'GENOMES/' ) )
gff_files = np.sort( gff_files[ np.array( [ f.endswith( '_genes.gff' ) for f in gff_files ] ) ] )
fna_files = np.array( os.listdir( 'GENOMES/' ) )
fna_files = np.sort( fna_files[ np.array( [ f.endswith( '.fna' ) for f in fna_files ] ) ] )

## temporary directory where intermediate files will be stored
##temp_dir = '/tmp/dvh-output/'   ## './dvh-output/'
temp_dir = tempfile.mktemp( dir='./' ) + '/'

sample_files = np.array( os.listdir( remote_dir ) )
sample_files = np.sort( sample_files[ np.array( [ f.startswith( 'Sample_' ) for f in sample_files ] ) ] )
if sample is not None:
    sample_files = sample_files[ np.where( sample_files == sample ) ]
print sample_files

if sample is None and len(sample_files) > 1:
    ut.system( "rm -rf dvh-output" )

try:
    os.makedirs( temp_dir )
    os.makedirs( './dvh-starIndex' )
    os.makedirs( './dvh-output' )
    #os.makedirs( remote_dir + 'dvh-starIndex' )
    #os.makedirs( remote_dir + 'dvh-output' )
except:
    pass

## note the genomeSAindexNbases option for smaller genomes, see:
## https://groups.google.com/d/msg/rna-star/j8KomjbDfW0/jsyL0c5m4xMJ
## also see here about the outSAMstrandField: http://seqanswers.com/forums/archive/index.php/t-30905.html
## update, do different index for each genome.
for org in genome_lookup.keys():
    out_dir = './dvh-starIndex/' + org ##fna.split('.')[0]
    print out_dir
    if os.path.exists( out_dir ):
        continue
    try:
        os.makedirs( out_dir )
    except:
        pass
    genomes = ' '.join( ['GENOMES/' + i + '.fna' for i in genome_lookup[org]] )
    cmd = ("STAR --runMode genomeGenerate --genomeSAindexNbases 7 --runThreadN 8 --genomeDir " + out_dir +
           " --genomeFastaFiles " + genomes)
    tmp = ut.system( cmd )

for SAMPLE in sample_files:
    #SAMPLE="Sample_US-1505885"

    print SAMPLE
    organisms = sample_tab.ix[sample_tab['Sample name']==SAMPLE].organisms.values[0]
    print 'ORGANISMS = ', organisms
    species = org_lookup[organisms]
    print 'SPECIES = ', species
    if len(species) <= 0:
        raise 'CANNOT FIND SPECIES FOR SAMPLE: ' + SAMPLE

    FASTQ1 = temp_dir + SAMPLE + "/" + SAMPLE + "_ALL1.fastq"
    FASTQ2 = temp_dir + SAMPLE + "/" + SAMPLE + "_ALL2.fastq"

    if os.path.isdir( remote_dir + 'dvh-output/' + SAMPLE ):
        #if os.path.isfile( remote_dir + 'dvh-output/' + SAMPLE + '/Aligned_sorted.out.bam' ):
        #if os.stat( remote_dir + 'dvh-output/' + SAMPLE + '/Aligned_sorted.out.bam' ).st_size > 0:
        continue

    try:
        #ut.system( "mkdir -p ./" + SAMPLE )
        os.makedirs( temp_dir + SAMPLE )
        os.makedirs( remote_dir + "/dvh-output/" + SAMPLE )
    except:
        pass

    FASTQ1 = FASTQ1 + ".gz"
    FASTQ2 = FASTQ2 + ".gz"
    ut.system( "rm -vf " + FASTQ1 + " " + FASTQ1 )

    fastq_files1 = glob.glob( remote_dir + SAMPLE + '/*_R1_001.fastq.gz' )
    fastq_files2 = glob.glob( remote_dir + SAMPLE + '/*_R2_001.fastq.gz' )

    ## Apparently you can just concatenate gzipped files together.... let's try it
    for i in range( len( fastq_files1 ) ):
        append = ' > ' if i == 0 else ' >> '
        ut.system( "cat " + fastq_files1[i] + append + FASTQ1 )
        ut.system( "cat " + fastq_files2[i] + append + FASTQ2 )

    if False:  ## check the reads using fastqc
        print ut.system( "gunzip -c " + FASTQ1 + "|wc -l" )
        print ut.system( "gunzip -c " + FASTQ2 + "|wc -l" )
        ut.system( "zcat " + FASTQ1 + " | fastqc -t 8 stdin" )
        #ut.system( "open stdin_fastqc.html" )

    ## see  here for outFilterMismatchNoverLmax option:
    ## http://seqanswers.com/forums/showthread.php?t=27645
    for org in species: ##genome_lookup.keys():
        print 'RUNNING: ', org
        try:
            os.makedirs( temp_dir + SAMPLE + "/" + org + "/" )
        except:
            pass

        genomeDir = "./dvh-starIndex/" + org + "/"
        outDir_temp = temp_dir + SAMPLE + "/" + org + "/"

        cmd = ("STAR --genomeDir " + genomeDir + " --runThreadN 7 --outSAMstrandField intronMotif "
           "--readFilesIn " +FASTQ1+ " " +FASTQ2+ " --readFilesCommand zcat "
               "--outFileNamePrefix " + outDir_temp)
        ##--outFilterMismatchNoverLmax 0.01 --outFilterMultimapNmax 100 --outFilterMultimapScoreRange 2 --outSAMstrandField None --outSAMmode Full --outSAMattributes Standard --outSAMunmapped None --outFilterType BySJout --outStd SAM
        ut.system( cmd )

        sam_file = outDir_temp + "Aligned.out.sam"
        ut.system( "pigz -fv9 " + sam_file )
        sam_filegz = outDir_temp + "Aligned.out.sam.gz"

        ## to read in R:
        ## x=fread(gzfile('dvh-output/Sample_US-1505885/Aligned.out.sam.gz'),skip=5)

        ## makes bam file and sorts; let it use more RAM than the default 768M
        ## see http://samtools.sourceforge.net/samtools.shtml
        ## also http://davetang.org/wiki/tiki-index.php?page=SAMTools#Converting_SAM_directly_to_a_sorted_BAM_file
        ## also http://biobits.org/samtools_primer.html
        ## Note -n is to sort by name instead of position, gets around the htseq-count buffer overflow error:
        ##   http://seqanswers.com/forums/showthread.php?t=41531
        ## Note I also added the bundle=True option to HTSeq.pair_SAM_alignments( read_seq, bundle=True ) in counts.py 
        ##    - this *should* count multiple alignments instead of throwing them away...?
        ## this doesn't work so I removed it again, for more info, look at 
        ##    http://www-huber.embl.de/users/anders/HTSeq/doc/counting.html#dealing-with-multiple-alignments 

        ## Now we're doing it for "stringtie", see here: 
        ## http://ccb.jhu.edu/software/stringtie/

        outDir_remote = remote_dir + "dvh-output/" + SAMPLE + "/" + org + "/"
        try:
            os.makedirs( outDir_remote )
        except:
            pass

        cmd = ("samtools view -bu -@ 8 " + sam_filegz + " | "
               "samtools sort -@ 8 -l 9 -m 2G - " + outDir_remote + "/Aligned_sorted")
        ut.system( cmd )
        bam_file = outDir_remote + "Aligned_sorted.bam"
        
        ##ut.system( "rm -fv dvh-output/" + SAMPLE + "/Aligned.out.sam.gz" )
        ##ut.system( "mv -fv " + temp_dir + SAMPLE + "/Aligned.out.sam.gz " + remote_dir + "dvh-output/" + SAMPLE )
        ##ut.system( "cp -fv " + temp_dir + SAMPLE + "/" + org + "/Aligned_sorted.bam " + 
        ##           remote_dir + "dvh-output/" + SAMPLE + "/" + org + "/" )

        ## to view bam file
        # samtools index dvh-output/${SAMPLE}/Aligned_sorted.out.bam
        # samtools tview dvh-output/${SAMPLE}/Aligned_sorted.out.bam

        ## I am removing all hits to multiple alignments which have NH:i:X, X>1, in this case most reads - 
        ##   this should speed up htseq-count which is really slow and htseq-count skips those anyway
        ## Update change cutoff to 2 - in case a read maps an orthologous gene e.g. in both DvH and G20 
        ## But need to replace NH:i:2 to NH:i:1 or it will still be ignored by htseq-count
        ## Then need to decrease min-aQual setting from default of 10. Jeez.
        ## Note quality score is int(-10.0*math.log10(1.0-1.0/Nmap)) so for 2 mappings, it's 3.
        ##"awk '($12==\"NH:i:1\"||$12==\"NH:i:2\"){print}' | "
        gffFiles = [ 'GENOMES/' + i + '_genes.gff' for i in genome_lookup[ org ] ]
        for gff in gffFiles:
            gff_b = os.path.basename(gff).replace( '_genes.gff', '' )
            out_file = outDir_remote + gff_b + '.counts'
            cmd = ("samtools view -@ 8 " + bam_file + " | "
                   "grep -E 'NH:i:[12]' | "
                   "sed 's/NH:i:2/NH:i:1/' | "
                   "htseq-count --mode=union --order=name --stranded=no --type=gene --idattr=locus_tag --minaqual=2 - "
                   + gff + " > " + out_file)
            ut.system( cmd )

        ## instead of htseq-count now lets use stringtie.
        ## change default args - shortest transcript 100 (instead of 200); max coverage 100,000,000 (instead of 1,000,000)
        ## can't do it, uses too much RAM
        # gffFiles = [ 'GENOMES/' + i + '_genes.gff' for i in genome_lookup[ org ] ]
        # for gff in gffFiles:
        #     gff_b = os.path.basename(gff).replace( '_genes.gff', '' )
        #     out_file = outDir_remote + gff_b + '.counts'
        #     cmd = ('stringtie ' + bam_file + ' -v -p 1 -m 200 -s 1000000 -e -G ' + gff + ' -o ' + out_file)
        #     ut.system( cmd )

        print "DONE:", SAMPLE, org
        ##ut.system( "rm -rvf " + remote_dir + "dvh-output/" + SAMPLE )
        ut.system( "cp -avf " + temp_dir + SAMPLE + "/" + org + " " + remote_dir + "dvh-output/" + SAMPLE )
        ##ut.system( "rm -rvf ./" + SAMPLE + "/" + org )    
        ut.system( "rm -rvf " + temp_dir + SAMPLE + "/" + org )

    ut.system( "rm -rvf " + temp_dir + SAMPLE )

## the --format=bam option requires pysam: pip install --user pysam
## but should be faster???. https://code.google.com/p/pysam/
#htseq-count --format=bam --mode=union --order=name --stranded=no --type=gene --idattr=locus_tag dvh-output/${SAMPLE}/Aligned_sorted.out.bam DvH_genes.gff > dvh-output/${SAMPLE}/DvH_ALL.counts

## or this. I had to edit HTSeq-0.6.1/HTSeq/scripts/count.py to add max_buffer_size=40000000 (40 million) parameter to HTSeq.pair_SAM_alignments_with_buffer() -- default is 3000000 (3 million)
#samtools view -bu -@ 8 dvh-output/${SAMPLE}/Aligned.out.sam.gz | samtools sort -@ 8 -l 9 -m 2G - dvh-output/${SAMPLE}/Aligned_sorted2.out
#samtools index dvh-output/${SAMPLE}/Aligned_sorted2.out.bam
#htseq-count --format=bam --mode=union --order=pos --stranded=no --type=gene --idattr=locus_tag dvh-output/${SAMPLE}/Aligned_sorted.out.bam DvH_genes.gff > dvh-output/${SAMPLE}/DvH_ALL.counts2

## or:
#cufflinks -p 8 -G DvH_genes.gff -u dvh-output/${SAMPLE}/Aligned_sorted2.out.bam -o dvh-output/${SAMPLE}/cufflinks_out

###################3###################3###################3###################3###################3
## TO USE BOWTIE INSTEAD OF STAR:
#bowtie2-build NC_002937.fna,NC_005863.fna dvh-bt2/dvh

#mkdir -p dvh-output/${SAMPLE}
#bowtie2 --very-sensitive -p 8 --reorder --mm  -1 ${FASTQ1} -2 ${FASTQ2} -x dvh-bt2/dvh -S dvh-output/${SAMPLE}/bt2_alignments.sam

#samtools view -bu -@ 8 dvh-output/${SAMPLE}/bt2_alignments.sam | samtools sort -@ 8 -l 9 -m 2000000000 - dvh-output/${SAMPLE}/bt2_alignments_sorted

#cufflinks -p 8 -G NC_002937_genes.gff -u dvh-output/${SAMPLE}/bt2_alignments_sorted.bam -o dvh-output/bt2_${SAMPLE}/
## END BOWTIE
