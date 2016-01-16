import sys,numpy

# 1. reading the data
print 'reading the data...'
# 1.1. read the FPKM matrix
inputFile='/Volumes/omics4tb/alomana/projects/dtp/results/genes.fpkm_table.v2.txt'
M=[]
geneNames=[]
with open(inputFile,'r') as f:
    header=f.readline()
    prelabels=header.split('\t')[1:]
    labels=[element.split('_')[0] for element in prelabels]
    f.next()
    for line in f:
        vector=line.split('\t')
        geneName=vector[0]
        geneNames.append(geneName)
        prevalues=vector[1:]
        values=[float(element) for element in prevalues]
        M.append(values)

# 1.2. reading the metadata
metaData={}
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.tsv'
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

print metaData

# 2. select the two controls and generate two output matrices
sampleControls300=[]; sampleControls1000=[]
for label in labels:
    if metaData[label]['epoch'] == 0 and metaData[label]['growth'] == 'exp' and metaData[label]['light'] == 'PM' and metaData[label]['co2'] == 300:
        sampleControls300.append(label)
    if metaData[label]['epoch'] == 0 and metaData[label]['growth'] == 'exp' and metaData[label]['light'] == 'PM' and metaData[label]['co2'] == 1000:
        sampleControls1000.append(label)

# 3. transforming the data 
print 'transforming the data...'
D=numpy.array(M)
F=D+1.

# 3.1. generating the 300 ppm relative expression matrix
selectedColumns=[labels.index(element) for element in sampleControls300]
C=F[:,selectedColumns]
Cmean=C.mean(axis=1)
R=F/Cmean[:,numpy.newaxis]
log2R=numpy.log2(R)

print 'saving the data for 300 ppm...'
outputFile='/Volumes/omics4tb/alomana/projects/dtp/results/relative.log2values.300.txt'
with open(outputFile,'w') as g:
    
    # writing the header
    g.write('tracking_id\t')
    for j in range(len(labels)):
        if metaData[labels[j]]['co2'] == 300 and labels[j] not in sampleControls300:
            g.write('%s\t'%labels[j])
    g.write('\n')

    # writing each line
    for i in range(len(geneNames)):
        g.write('%s\t'%geneNames[i])
        for j in range(len(labels)):
            if metaData[labels[j]]['co2'] == 300 and labels[j] not in sampleControls300:
                value=log2R[i,labels.index(labels[j])]
                g.write('%s\t'%value)
        g.write('\n')

## print 'F'
## print F[120:130,:]
## print 'Cmean'
## print Cmean[120:130]
## print 'R'
## print R[120:130,:]
## print 'log2R'
## print log2R[120:130,:]

# 3.2. generating the 1,000 ppm relative expression matrix
selectedColumns=[labels.index(element) for element in sampleControls1000]
C=F[:,selectedColumns]
Cmean=C.mean(axis=1)
R=F/Cmean[:,numpy.newaxis]
log2R=numpy.log2(R)

print 'saving the data for 1,000 ppm...'
outputFile='/Volumes/omics4tb/alomana/projects/dtp/results/relative.log2values.1000.txt'
with open(outputFile,'w') as g:
    
    # writing the header
    g.write('tracking_id\t')
    for j in range(len(labels)):
        if metaData[labels[j]]['co2'] == 1000 and labels[j] not in sampleControls1000:
            g.write('%s\t'%labels[j])
    g.write('\n')

    # writing each line
    for i in range(len(geneNames)):
        g.write('%s\t'%geneNames[i])
        for j in range(len(labels)):
            if metaData[labels[j]]['co2'] == 1000 and labels[j] not in sampleControls1000:
                value=log2R[i,labels.index(labels[j])]
                g.write('%s\t'%value)
        g.write('\n')

# 4. final message
print '... done.'
