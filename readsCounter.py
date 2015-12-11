import sys,os

'''
this script quantifies the read counts using HTSeq. To be run in osiris.
'''

def htseqCounter(folder):

    '''
    this function calls HTSeq
    '''

    htseqExecutable='python -m HTSeq.scripts.count -r pos -f bam -t gene -s reverse '
    bamFile=bamFilesDir+folder+'/Aligned.sortedByCoord.out.bam '
    genomeAnnotationFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
    outputDirection=' > %s%s.txt'%(countsDir,folder)
    
    cmd=htseqExecutable+bamFile+genomeAnnotationFile+outputDirection

    print
    print cmd
    print
    os.system(cmd)

    #sys.exit()

    return None

# 0. defining user variables
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/reads/tippingPoints/bamFiles/'
countsDir='/proj/omics4tb/alomana/projects/dtp/data/reads/tippingPoints/counts/'

# 1. defining the BAM files
bamRoots=os.listdir(bamFilesDir)

# 2. calling HTSeq
for folder in bamRoots:
    htseqCounter(folder)
