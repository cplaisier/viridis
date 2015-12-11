import os,sys

'''
this script finds the clean FASTQ files and calls STAR for the reads alignment.
'''

def genomeIndexer():

    '''
    this function creates the genome index. Should be run only once.
    '''

    flag1=' --runMode genomeGenerate'
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --genomeDir %s'%genomeIndexDir
    flag4=' --genomeFastaFiles %s'%genomeFastaFile
    flag5=' --sjdbGTFfile %s'%genomeAnnotationFile
    flag6=' --sjdbGTFtagExonParentTranscript Parent --sjdbOverhang 150'

    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag5+flag6

    print
    print cmd
    print
    os.system(cmd)

    return None

def STARcalling(inputFile):

    '''
    this function calls STAR
    '''
    
    finalDir=bamFilesDir+inputFile+'/'
    if os.path.exists(finalDir) == False:
        os.mkdir(finalDir)

    readF1=readsFilesDir+inputFile+'_L001_R1_001.clean.fastq.gz'
    readF2=readsFilesDir+inputFile+'_L002_R1_001.clean.fastq.gz'
    readF3=readsFilesDir+inputFile+'_L003_R1_001.clean.fastq.gz'
    readF4=readsFilesDir+inputFile+'_L004_R1_001.clean.fastq.gz'

    readR1=readsFilesDir+inputFile+'_L001_R2_001.clean.fastq.gz'
    readR2=readsFilesDir+inputFile+'_L002_R2_001.clean.fastq.gz'
    readR3=readsFilesDir+inputFile+'_L003_R2_001.clean.fastq.gz'
    readR4=readsFilesDir+inputFile+'_L004_R2_001.clean.fastq.gz'

    flag1=' --genomeDir %s'%genomeIndexDir
    flag2=' --runThreadN %s'%numberOfThreads
    flag3=' --readFilesIn %s,%s,%s,%s %s,%s,%s,%s'%(readF1,readF2,readF3,readF4,readR1,readR2,readR3,readR4)    
    flag4=' --outFileNamePrefix %s'%finalDir
    flag5=' --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c --quantMode TranscriptomeSAM GeneCounts'

    cmd=STARexecutable+flag1+flag2+flag3+flag4+flag5
    
    print
    print cmd
    print
    os.system(cmd)
    
    return None

# 0. defining several input/output paths
readsFilesDir='/proj/omics4tb/alomana/projects/dtp/data/reads/tippingPoints/cleanReads/'
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/reads/tippingPoints/bamFiles/'
STARexecutable='/proj/omics4tb/alomana/software/STAR-master/bin/Linux_x86_64/STAR'
genomeIndexDir='/proj/omics4tb/alomana/projects/dtp/data/genomeIndex'
genomeFastaFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.dna.genome.fa'
genomeAnnotationFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
numberOfThreads=8

# 1. recover the clean FASTQ files
print 'reading FASTQ files...'
allTags=[]
allFiles=os.listdir(readsFilesDir)
for element in allFiles:
    tag=element.split('_L00')[0]
    allTags.append(tag)
inputFiles=list(set(allTags))

# 2. making genome indexes
#print 'making genome index...'
#genomeIndexer()
       
# 3. calling STAR
print 'calling STAR...'
for inputFile in inputFiles:
    STARcalling(inputFile)
