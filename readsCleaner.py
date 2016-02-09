### this script calls Trimmomatic to clean the reads
import os,sys

def trimmomaticCaller(instance):
    '''
    This function deals with the trimming of the reads using Trimmomatic. Recommended options, ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    '''

    print 'working with file',instance

    logFile=logFilesDir+instance+'.messagesFromTrimming.txt'

    inputFile1=rawReadsDir+instance+'_R1_001.fastq.gz'
    inputFile2=rawReadsDir+instance+'_R2_001.fastq.gz'

    outputFile1a=cleanReadsDir+instance+'_R1_001.clean.fastq.gz'
    outputFile1b=cleanReadsDir+instance+'_R1_001.garbage.fastq.gz'

    outputFile2a=cleanReadsDir+instance+'_R2_001.clean.fastq.gz'
    outputFile2b=cleanReadsDir+instance+'_R2_001.garbage.fastq.gz'

    cmd='/Library/Internet\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /Users/alomana/software/Trimmomatic-0.33/trimmomatic-0.33.jar PE -threads 4 -phred33 -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(logFile,inputFile1,inputFile2,outputFile1a,outputFile1b,outputFile2a,outputFile2b,path2Adapter)  

    print cmd
    os.system(cmd)
    print
    
    return None

# 0. defining user variables
rawReadsDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/rawData/secondRun/'
cleanReadsDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cleanReads/secondRun/'
logFilesDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/logFilesTrimmomatic/'

path2Adapter='/Users/alomana/projects/ap/seqs/src/adapters/TruSeq3-PE-2.fa'

# 1. reading the files
foundFiles=os.listdir(rawReadsDir)
readsFiles=[element for element in foundFiles if not element.startswith('.')]

allFiles=[]
for readsFile in readsFiles:
    label=readsFile.split('/')[-1].split('_R')[0]
    allFiles.append(label)
sequencingInstances=list(set(allFiles))

# 2. calling Trimmomatic
for instance in sequencingInstances:
    trimmomaticCaller(instance)
