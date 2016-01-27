### this script locates epoch 2 samples into a light vs. growth space
import sys,numpy,matplotlib
from matplotlib import pyplot

def boxPlotGrapher(classifiers,flag):

    '''
    this function plots the expression of the classifiers in log10 scale using boxplots
    '''

    # 1. defining the expression values
    if flag == 'light':
        expressionA=expressionRetriever(classifiers,flag,'AM')
        expressionB=expressionRetriever(classifiers,flag,'PM')
    elif flag == 'growth':
        expressionA=expressionRetriever(classifiers,flag,'exp')
        expressionB=expressionRetriever(classifiers,flag,'sta')
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # 3. plotting
    boxPlotPosition=0
    listOfClassifiers=sorted(classifiers,key=classifiers.__getitem__,reverse=True)# ranking the classifiers
    for geneID in listOfClassifiers: 
        boxPlotPosition=boxPlotPosition+1
        # 3.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        # 3.2. actual plotting
        bp=matplotlib.pyplot.boxplot([logx],positions=[boxPlotPosition],patch_artist=True)
        setBoxColors(bp,'orange')
        bp=matplotlib.pyplot.boxplot([logy],positions=[boxPlotPosition],patch_artist=True)
        setBoxColors(bp,'green')
        
    # 3.3. closing the figure
    matplotlib.pyplot.xlim([0,boxPlotPosition+1])
    matplotlib.pyplot.ylim([-0.2,3.75])
    theXticks=range(boxPlotPosition)
    theXticksPosition=[element+1 for element in theXticks]
    matplotlib.pyplot.xticks(theXticksPosition,listOfClassifiers,rotation=-90)
    matplotlib.pyplot.ylabel('log10 FPKM')
    matplotlib.pyplot.tight_layout(pad=0.5)
    matplotlib.pyplot.savefig('boxplots_%s.pdf'%flag)
    matplotlib.pyplot.clf()

    return None

def classifiersFilter(classifiers,flag):

    '''
    this function selects the classifiers that have at least 0.1 separation between max and min values into the two cases in log space
    '''
    
    # 1. defining the expression values
    if flag == 'light':
        selectedExpressionA=expressionRetriever(classifiers,flag,'AM')
        selectedExpressionB=expressionRetriever(classifiers,flag,'PM')
    elif flag == 'growth':
        selectedExpressionA=expressionRetriever(classifiers,flag,'exp')
        selectedExpressionB=expressionRetriever(classifiers,flag,'sta')
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # 2. computing the distances
    separations={}
    for geneID in classifiers: 
        
        # 2.1 transforming into log10 scale
        x=numpy.array(selectedExpressionA[geneID])
        y=numpy.array(selectedExpressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        # defining the separation
        if numpy.mean(logx) > numpy.mean(logy):
            separation=numpy.min(logx)-numpy.max(logy)
        elif numpy.mean(logy) > numpy.mean(logx):
            separation=numpy.min(logy)-numpy.max(logx)
        else:
            print 'error computing the separation between the two distributions from classifiersFilter'
            sys.exit()
        if separation > 0.:
            separations[geneID]=separation

    return separations

def classifiersRetriever(flag):

    '''
    this function reads the output of cuffdiff and select the genes that pass the following rule: 1. log2 fold > 1., and 2. statistical significance
    '''

    classifiersFolds={}
    classifiersSignificances={}
    
    inputFile=cuffdiffDir+'%s/gene_exp.diff'%flag
    with open(inputFile,'r') as f:
        f.next()
        for line in f:
            vector=line.split('\t')
            geneName=vector[0]
            significance=vector[-1].replace('\n','')
            foldChange=abs(float(vector[-5]))
            qValue=float(vector[-2])
            if significance == 'yes' and foldChange > 1.: # the selection rule is log2 fold > 1.0 and statistically different
                classifiersFolds[geneName]=foldChange
                classifiersSignificances[geneName]=qValue
    # sorting classifiers
    listOfClassifiers=sorted(classifiersFolds,key=classifiersFolds.__getitem__,reverse=True)

    return listOfClassifiers

def expressionReader():

    '''
    this function reads the matrix of expression in FPKM into the format of a dictionary
    '''

    expression={}
    with open(expressionFile,'r') as f:
        header=f.readline()
        prelabels=header.split('\t')[1:]
        labels=[element.split('_')[0] for element in prelabels]
        f.next()
        for line in f:
            vector=line.split('\t')

            geneName=vector[0]
            expression[geneName]={}
            
            preValues=vector[1:]
            values=[float(element) for element in preValues]
            for i in range(len(values)):
                expression[geneName][labels[i]]=values[i]

    return expression

def expressionRetriever(classifiers,flag,condition):

    '''
    this function returns the expression values for a set of genes under a specific condition
    '''

    # 1. selecting the acceptable samples
    selectedSamples=[]
    for sampleID in metaData.keys():
        if metaData[sampleID]['epoch'] != 2:
            if metaData[sampleID][flag] == condition:
                selectedSamples.append(sampleID)
                
    # 2. defining the expression values
    selectedExpression={}
    for geneName in classifiers:
        selectedExpression[geneName]=[]
        for sampleID in selectedSamples:
            value=expression[geneName][sampleID]
            selectedExpression[geneName].append(value)

    return selectedExpression

def metadataReader():

    '''
    this function returns a dictionary with the metadata of the expression data
    '''
    
    metaData={}

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

    return metaData

def newSpaceMapper(sampleID):

    '''
    this function returns a sample into the new space coordinates
    '''

    x=numpy.random.random()
    y=numpy.random.random()

    return x,y

def setBoxColors(bp,theColor):

    '''
    this function access the elements of a boxplot and colors them appropriately
    '''

    matplotlib.pyplot.setp(bp['boxes'],color=theColor)
    matplotlib.pyplot.setp(bp['caps'],color=theColor)
    matplotlib.pyplot.setp(bp['whiskers'],color=theColor,ls='-')
    matplotlib.pyplot.setp(bp['fliers'],markeredgecolor=theColor,marker='+')
    matplotlib.pyplot.setp(bp['medians'],color=theColor)    

    return None

### MAIN

# 0. preliminaries
print 'initializing variables...'
# 0.1. user defined variables and paths
cuffdiffDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.tsv'

# 0.2. reading metadata
metaData=metadataReader()

# 0.3. reading expression
expression=expressionReader()

# 1. recover the DET genes. Rank them on q-value and plot the boxplots. select the best ones.
print 'recovering classifiers...'

lightClassifiers=classifiersRetriever('light')
print len(lightClassifiers),'light classifiers detected.'

growthClassifiers=classifiersRetriever('growth')
print len(growthClassifiers),'growth classifiers detected.'

# 1.2. filtering the classifiers based on separation
print 'filtering classifiers based on separation...'
lightFilteredClassifiers=classifiersFilter(lightClassifiers,'light')
growthFilteredClassifiers=classifiersFilter(growthClassifiers,'growth')
print len(lightFilteredClassifiers),'filtered light classifiers.'
print len(growthFilteredClassifiers),'filtered growth classifiers.'

# 1.3. plotting a boxplots of the best classifiers
print 'plotting expression graphs for classifiers...'
boxPlotGrapher(lightFilteredClassifiers,'light')
boxPlotGrapher(growthFilteredClassifiers,'growth')

# 2. map samples into a new space of dark/light distributed in x:-2:-1/1:2 and stationary/exponential y:-2:-1/1:2
print 'mapping samples into new space...'
for sampleID in metaData.keys():
    x,y=newSpaceMapper(sampleID)
    matplotlib.pyplot.plot(x,y,'ok')
    matplotlib.pyplot.tight_layout(pad=0.5)
    matplotlib.pyplot.savefig('sampleLocation.pdf')
    matplotlib.pyplot.clf()
    #if metaData[sampleID]['epoch'] == 2:

# 3. final message
print '... analysis completed.'
        
