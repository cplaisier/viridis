### this script finds the DETs specifically for light/dark and exponential/stationary phases 
import os,sys

# 0. defining some variables
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
cufflinksDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/'
maskFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/mask.gff3'
outputDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
numberOfThreads=12

# 1. reading the metadata
metaData={}
metaDataFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.tsv'
with open(metaDataFile,'r') as f:
    f.next()
    for line in f:
        vector=line.split('\t')
        if vector[4] != '':
            sampleID=vector[5]
            if vector[0] != '':
                epoch=int(vector[0])
            growth=vector[1]
            light=vector[2]
            co2=int(vector[3].replace(',',''))
            replicate=vector[4]
            metaData[sampleID]={}
            metaData[sampleID]['growth']=growth
            metaData[sampleID]['epoch']=epoch
            metaData[sampleID]['light']=light
            metaData[sampleID]['co2']=co2
            metaData[sampleID]['replicate']=replicate

# 2. select the samples. excluding epoch 2 of 1,000 ppm
lightSamples=[]
darkSamples=[]
exponentialSamples=[]
stationarySamples=[]

for sampleID in metaData:
    if metaData[sampleID]['epoch'] != 2:
        # 1.1. select light/dark samples
        if metaData[sampleID]['light'] == 'AM':
            lightSamples.append(sampleID)
        elif metaData[sampleID]['light'] == 'PM':
            darkSamples.append(sampleID)
        else:
            print 'error at light samples selection'
            sys.exit()

        # 1.2. select exponential/stationary samples
        if metaData[sampleID]['growth'] == 'exp':
            exponentialSamples.append(sampleID)
        elif metaData[sampleID]['growth'] == 'sta':
            stationarySamples.append(sampleID)
        else:
            print 'error at growth samples selection'
            sys.exit()

# 3. run cuffdiff
term1='cuffdiff %s '%(gtfFile)
term3='-p %s '%numberOfThreads
term4='-M %s '%maskFile
term5='--library-type fr-firststrand '
term6='--multi-read-correct '

paths=os.listdir(cufflinksDir)
cxbPaths=[cufflinksDir+path+'/abundances.cxb' for path in paths]

# 3.1 running cuffdiff for light/dark samples
term2='-o %s/light '%outputDir

pathsA=[]
pathsB=[]
for sampleID in lightSamples:
    for folder in cxbPaths:
        if sampleID in folder:
            pathsA.append(folder)
for sampleID in darkSamples:
    for folder in cxbPaths:
        if sampleID in folder:
            pathsB.append(folder)
term7='%s %s'%(','.join(pathsA),','.join(pathsB))

cmd=term1+term2+term3+term4+term5+term6+term7

print
print cmd
print

os.system(cmd)

# 3.2 running cuffdiff for exp/sta samples
term2='-o %s/growth '%outputDir

pathsA=[]
pathsB=[]
for sampleID in exponentialSamples:
    for folder in cxbPaths:
        if sampleID in folder:
            pathsA.append(folder)
for sampleID in stationarySamples:
    for folder in cxbPaths:
        if sampleID in folder:
            pathsB.append(folder)
print len(pathsA),len(pathsB)
term7='%s %s'%(','.join(pathsA),','.join(pathsB))

cmd=term1+term2+term3+term4+term5+term6+term7

print
print cmd
print

os.system(cmd)
