import sys,os

def cuffnormCaller():

    '''
    this function calls cuffnorm using all generated abundance binary files by cuffquant
    '''

    outputDir=cufflinksDir+'allSamples'

    term1='cuffnorm %s '%(gtfFile)
    term2=' '.join(abundanceFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    term6='-L '+','.join(labels)+' '
    cmd=term1+term2+term3+term4+term5+term6
    
    print
    print cmd
    print

    os.system(cmd)

    return None

def cuffquantCaller(inputFile):

    '''
    this function calls the different steps of the cufflinks pipeline.
    '''

    label=inputFile.split('/')[-2]
    outputDir=cufflinksDir+label
    
    # cuffquant
    term1='cuffquant %s '%(gtfFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='-M %s '%maskFile
    term5='--library-type fr-firststrand '
    term6='--multi-read-correct '
    cmd=term1+term2+term3+term4+term5+term6+inputFile

    print
    print cmd
    print

    os.system(cmd)

    return None

# 0. defining input files
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cufflinksDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/'
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
maskFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/mask.gff3'
numberOfThreads=12

# 1. defining the BAM and abundance files
roots=os.listdir(bamFilesDir)

###
roots.remove('secondRun')
###

bamFiles=[bamFilesDir+element+'/Aligned.sortedByCoord.out.bam' for element in roots]
abundanceFiles=[cufflinksDir+element+'/abundances.cxb' for element in roots]
labels=[element.split('_')[-1] for element in roots]

# 2. calling cuffquantCaller 
print 'calling cuffquant...'
for inputFile in bamFiles:
    cuffquantCaller(inputFile)

# 3. calling cuffnorm
print 'calling cuffnorm...'
cuffnormCaller()
