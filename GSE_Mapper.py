### this script locates all epochs samples into a light vs. growth space characterized in epoch=0

import sys,numpy,matplotlib,random
from matplotlib import pyplot
from matplotlib.patches import Ellipse
from previousDescriptors import *

def boxPlotGrapher(classifiers,borders,flag):

    '''
    this function plots the expression of the classifiers in log10 scale using boxplots
    '''

    # 0. defining a variable for recapitulating the data borders necessary for mapping samples into new space
    borders[flag]={}
    
    # 1. defining the expression values
    if flag == 'diurnal':
        expressionA=descriptorsExpressionRetriever(classifiers,flag,'light')
        expressionB=descriptorsExpressionRetriever(classifiers,flag,'dark')           
    elif flag == 'growth':
        expressionA=descriptorsExpressionRetriever(classifiers,flag,[1,2])
        expressionB=descriptorsExpressionRetriever(classifiers,flag,[4,5])
    else:
        print flag
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # 3. plotting and recovering the info for the mapping samples
    boxPlotPosition=0
    listOfClassifiers=sorted(classifiers,key=classifiers.__getitem__,reverse=True) # ranking the classifiers
    for geneID in listOfClassifiers:
        boxPlotPosition=boxPlotPosition+1
        borders[flag][geneID]={}
        
        # 3.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        # 3.2. actual plotting
        if boxplotPlotting == True:
            bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
            setBoxColors(bp,'orange')
            bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
            setBoxColors(bp,'darkgreen')

        # 3.3. saving the info for the mapping samples
        xa=numpy.min(x); xb=numpy.median(x); xc=numpy.max(x); sdx=numpy.std(x)
        ya=numpy.min(y); yb=numpy.median(y); yc=numpy.max(y); sdy=numpy.std(y)
        if xb > yb:
            center=((xa-yc)/2.)+yc
        else:
            center=((ya-xc)/2.)+xc
        # incorporating the data
        #! w=weightRankCalculator(geneID,listOfClassifiers)
        w=weightNMLCalculator(geneID,flag)
        if flag == 'diurnal':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        elif flag == 'growth':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        else:
            print 'error handling flags from boxPlotGrapher (bis)'
            sys.exit()
        
    # 3.3. closing the figure
    if boxplotPlotting == True:
        matplotlib.pyplot.xlim([0,boxPlotPosition+1])
        #matplotlib.pyplot.ylim([-0.2,5.])
        theXticks=range(boxPlotPosition)
        theXticksPosition=[element+1 for element in theXticks]
        matplotlib.pyplot.xticks(theXticksPosition,listOfClassifiers,rotation=-90,fontsize=2)
        matplotlib.pyplot.ylabel('log2 fold')
        matplotlib.pyplot.tight_layout(pad=2.5)
        #! matplotlib.pyplot.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.33)
        matplotlib.pyplot.savefig('boxplots_%s.pdf'%flag)
        matplotlib.pyplot.clf()

    return borders

def classifiersFilter(classifiers,flag):

    '''
    this function filters out descriptors
    '''

    # 1. defining the expression values
    if flag == 'diurnal':
        expressionA=descriptorsExpressionRetriever(classifiers,flag,'light')
        expressionB=descriptorsExpressionRetriever(classifiers,flag,'dark')        
    elif flag == 'growth':
        expressionA=descriptorsExpressionRetriever(classifiers,flag,[1,2])
        expressionB=descriptorsExpressionRetriever(classifiers,flag,[4,5])
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # 2. computing the separations
    separations={}
    for geneID in classifiers: 
        
        # 2.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        # defining the separation
        if numpy.mean(x) > numpy.mean(y):
            separation=numpy.min(x)-numpy.max(y)
        elif numpy.mean(y) > numpy.mean(x):
            separation=numpy.min(y)-numpy.max(x)
        else:
            print 'error computing the separation between the two distributions from classifiersFilter'
            sys.exit()
        if separation > 0.:
            separations[geneID]=separation

    return separations

def classifiersRetriever(flag):

    '''
    this function matches the previously found descriptors with the array data
    '''
    
    if flag == 'diurnal':
        fpkmDescriptors=diurnalFPKMdescriptors
    elif flag == 'growth':
        fpkmDescriptors=growthFPKMdescriptors

    listOfClassifiers=[]
    for descriptor in fpkmDescriptors:
        for putative in expression.keys():
            if 'draft' in descriptor:
                a=int(descriptor.split('gene:Thapsdraft')[1])
            else:
                a=int(descriptor.split('gene:Thaps')[1])
            label='THAPSDRAFT_'
            if label in putative:
                b=int(putative.split('THAPSDRAFT_')[1])
            else:
                b=''
            if a == b:
                listOfClassifiers.append(putative)

    return listOfClassifiers

def expressionReader():

    '''
    this function reads the matrix of expression in FPKM into the format of a dictionary
    '''

    expression={}
    with open(expressionFile,'r') as f:
        header=f.readline()
        prelabels=header.split('\t')
        labels=[element.split('_')[0] for element in prelabels[11:30]]
        for line in f:
            vector=line.split('\t')
            geneName=vector[1]
            expression[geneName]={}
            
            preValues=vector[11:30]
            values=[float(element) for element in preValues]
            for i in range(len(values)):
                expression[geneName][labels[i]]=values[i]

    return expression

def descriptorsExpressionRetriever(classifiers,flag,condition):

    '''
    this function returns the expression values for a set of genes under a specific condition
    '''

    # 1. selecting the acceptable samples
    selectedSamples=[]
    for sampleID in metaData.keys():
        if metaData[sampleID]['co2'] == 400:
            if condition == 'light' or condition == 'dark':
                if metaData[sampleID][flag] == condition:
                    selectedSamples.append(sampleID)
            if len(condition) == 2:
                if metaData[sampleID][flag] in condition:
                    selectedSamples.append(sampleID)

    print flag,condition,selectedSamples
    # 2. defining the expression values
    selectedExpression={}
    for geneName in classifiers:
        selectedExpression[geneName]=[]
        for sampleID in selectedSamples:
            value=expression[geneName][sampleID]
            selectedExpression[geneName].append(value)

    return selectedExpression

def loadCalculator(sampleID,flag):

    '''
    this function computes the value of the sample in the new axis 
    '''

    averageLoad=0.
    
    for classifier in borders[flag].keys():
        s=expression[classifier][sampleID] # sample expression level

        a=borders[flag][classifier][0] # min of positive
        b=borders[flag][classifier][1] # median of positive
        c=borders[flag][classifier][2] # max of positive

        d=borders[flag][classifier][3] # min of negative
        e=borders[flag][classifier][4] # median of negative
        f=borders[flag][classifier][5] # max of negative

        NML=borders[flag][classifier][6] # No Man's Land
        w=borders[flag][classifier][7] # weight
        
        # differentiate if it's positive (AM,exp) or negative (PM,sta) or miss-regulated
        positive=None
        missR=None
        if b > e:
            if s > NML:
                positive=True
            else:
                positive=False
            if s < a and s > f:
                missR=True
            else:
                missR=False
        if e > b:
            if s > NML:
                positive=False
            else:
                positive=True
            if s < d and s > c:
                missR=True
            else:
                missR=False
        
        # working with the misregulated samples
        if missR == True:
            if b > e:
                if s > NML:
                    stretch=abs(a-NML)
                    value=(s-NML)/stretch
                else:
                    stretch=abs(f-NML)
                    value=-abs(s-NML)/stretch
            else:
                if s > NML:
                    stretch=abs(d-NML)
                    value=-abs(s-NML)/stretch
                else:
                    stretch=abs(c-NML)
                    value=abs(s-NML)/stretch
                        
        # dealing with values within previously observed
        else:
            if positive == True: # it is a light sample
                # assuming light boxplot is above
                if b > e:
                    if s > b:
                        stretch=abs(c-b)
                        deviation=abs(s-b)
                        value=1.5+0.5*(deviation/stretch)
                    else:
                        stretch=abs(b-a)
                        deviation=abs(s-b)
                        value=1.5-0.5*(deviation/stretch)
                # assuming light boxplot is below
                else:
                    if s > b:
                        stretch=abs(c-b)
                        deviation=abs(s-b)
                        value=1.5-0.5*(deviation/stretch)
                    else:
                        stretch=abs(b-a)
                        deviation=abs(s-b)
                        if stretch != 0.:
                            value=1.5+0.5*(deviation/stretch)
                        else:
                            value=2.
            else: # it is a dark sample
                # assuming light boxtplot is above
                if b > e:
                    if s > e:
                        stretch=abs(f-e)
                        deviation=abs(s-e)
                        value=-1.5+0.5*(deviation/stretch)
                    else:
                        stretch=abs(e-d)
                        deviation=abs(s-e)
                        value=-1.5-0.5*(deviation/stretch)
                # assuming light boxplot is below
                else:
                    if s > e:
                        stretch=abs(f-e)
                        deviation=abs(s-e)
                        value=-1.5-0.5*(deviation/stretch)
                    else:
                        stretch=abs(e-d)
                        deviation=abs(s-e)
                        value=-1.5+0.5*(deviation/stretch)

        # weighting the value
        term=w*value
        averageLoad=averageLoad+term

    return averageLoad

def metadataReader():

    '''
    this function returns a dictionary with the metadata of the expression data
    '''
    
    metaData={}

    with open(expressionFile,'r') as f:
        header=f.readline()
        vector=header.split('\t')

        for sampleID in vector:
            if 'X' in sampleID:

                # co2 levels
                co2=int(sampleID.split('.')[0].replace('X',''))

                # diurnal
                if sampleID.split('.')[-1] == 'Dk':
                    diurnal='dark'
                elif sampleID.split('.')[-1] == 'Lt':
                    diurnal='light'

                # day
                growth=int(sampleID.split('.')[1].replace('Day',''))

                # filling up
                metaData[sampleID]={}
                metaData[sampleID]['co2']=co2
                metaData[sampleID]['diurnal']=diurnal
                metaData[sampleID]['growth']=growth

    return metaData

def newCoordinateCalculator(sampleID):

    '''
    this function computes the new coordinates of a sample given the borders calculated from the filtered classifiers
    '''

    x=-loadCalculator(sampleID,'diurnal')
    y=loadCalculator(sampleID,'growth')
    
    return x,y

def newSpaceMapper(flag):

    '''
    this function plots the samples into a new space
    '''

    preselectedSamples=[sampleID for sampleID in metaData.keys() if metaData[sampleID]['co2'] == int(flag)]
    print 'selected ', len(preselectedSamples), 'samples for plotting on ',flag, 'condition.'

    # starting the figure
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(111)

    for sampleID in preselectedSamples:
        
        x,y=newCoordinateCalculator(sampleID)
        print sampleID,x,y
        
        theSize=10
        theAlpha=.85
        theMarker='o'
        theColor='black'

        theDay=sampleID.split('.')[1].replace('Day','')

        ax.plot(x,y,marker=theMarker,mew=1,color=theColor,ms=theSize,alpha=theAlpha,mec='None',zorder=10)
        ax.plot(x,y+0.125,marker="$%s$"%theDay)

    # connecting the dots
    sortedSamples=['X800.Day1.Lt', 'X800.Day1.Dk', 'X800.Day2.Lt', 'X800.Day2.Dk', 'X800.Day3.Lt', 'X800.Day3.Dk', 'X800.Day4.Lt', 'X800.Day4.Dk', 'X800.Day5.Dk']
    
    for i in range(len(sortedSamples)-1):
        initial=sortedSamples[i]
        last=sortedSamples[i+1]
        a,b=newCoordinateCalculator(initial)
        c,d=newCoordinateCalculator(last)
        ax.plot([a,c],[b,d],':k',zorder=5,alpha=0.75)
    

    # setting ranges
    matplotlib.pyplot.xlim([-2.,2.])
    matplotlib.pyplot.ylim([-2.,2.])
    matplotlib.pyplot.xticks([-1.5,1.5],['light','dark'],fontsize=24)
    matplotlib.pyplot.yticks([-1.5,1.5],['late','early'],fontsize=24)
    matplotlib.pyplot.xlabel('diurnal phase',fontsize=24)
    matplotlib.pyplot.ylabel('growth phase',fontsize=24)
    matplotlib.pyplot.tight_layout(pad=2.5)

    # defining health zones
    matplotlib.pyplot.plot([-2,-1],[1,1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([-1,-1],[1,2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([2,1],[1,1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([1,1],[1,2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([2,1],[-1,-1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([1,1],[-1,-2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([-2,-1],[-1,-1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([-1,-1],[-1,-2],color='black',alpha=0.2)

    # defining the grid lines
    matplotlib.pyplot.plot([-1.5,-1.5],[1,2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([-2,-1],[1.5,1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([-1.5,-1.5],[-1,-2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([-2,-1],[-1.5,-1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([1.5,1.5],[-1,-2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([2,1],[-1.5,-1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([1.5,1.5],[1,2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([2,1],[1.5,1.5],color='black',ls=':',alpha=0.4)

    # defining the misregulation zone
    matplotlib.pyplot.plot([-1,-1],[-1,1],color='magenta',alpha=0.5,lw=2.)
    matplotlib.pyplot.plot([-1,1],[1,1],color='magenta',alpha=0.5,lw=2.)
    matplotlib.pyplot.plot([1,1],[1,-1],color='magenta',alpha=0.5,lw=2.)
    matplotlib.pyplot.plot([-1,1],[-1,-1],color='magenta',alpha=0.5,lw=2.)

    matplotlib.pyplot.savefig('sampleLocation.%s.pdf'%(str(int(flag))))
    matplotlib.pyplot.clf()
       
    return None

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

def weightNMLCalculator(geneID,flag):

    '''
    this function computes the weight of the classifier based on the empty space between the classifiers
    '''

    if flag == 'diurnal':
        sumSpaces=sum(diurnalFilteredClassifiers.values())
        value=diurnalFilteredClassifiers[geneID]
    elif flag == 'growth':
        sumSpaces=sum(growthFilteredClassifiers.values())
        value=growthFilteredClassifiers[geneID]
    else:
        print 'error from weightNMLCalculator'
        sys.exit()

    weight=value/sumSpaces

    return weight

### MAIN

# use the points defined in day 1 and use the gap as unit of measurement. rank equal weights.

# 0. preliminaries
print 'initializing variables...'

# 0.1. user defined variables and paths
expressionFile='../data/GSE45252_value_Matrix.txt'
boxplotPlotting=False

# 0.2. reading metadata
metaData=metadataReader()

# 0.3. reading expression
expression=expressionReader()

# 1. define descriptors

# 1.1. recover descriptors
print 'recovering descriptors...'

diurnalClassifiers=classifiersRetriever('diurnal')
growthClassifiers=classifiersRetriever('growth')

print len(diurnalClassifiers),'light descriptors found.'
print len(growthClassifiers),'growth descriptors found.'

# 1.2. filtering the classifiers based on separation
print 'filtering classifiers based on separation...'

diurnalFilteredClassifiers=classifiersFilter(diurnalClassifiers,'diurnal')
growthFilteredClassifiers=classifiersFilter(growthClassifiers,'growth')
print len(diurnalFilteredClassifiers),'filtered light classifiers.'
print len(growthFilteredClassifiers),'filtered growth classifiers.'

# 1.3. plotting a boxplots of the best classifiers
print 'plotting expression graphs for classifiers...'
borders={}
borders=boxPlotGrapher(diurnalFilteredClassifiers,borders,'diurnal')
borders=boxPlotGrapher(growthFilteredClassifiers,borders,'growth')

# 2. map samples into a new space of dark/light distributed in x:-2:-1/1:2 and stationary/exponential y:-2:-1/1:2
print 'mapping samples into new space...'
newSpaceMapper('800')

# 4. final message
print '... analysis completed.'
        