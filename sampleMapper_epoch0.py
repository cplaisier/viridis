### this script locates all epochs samples into a light vs. growth space characterized in epoch=0

import sys,numpy,matplotlib,random
from matplotlib import pyplot
from matplotlib.patches import Ellipse

def addArrows():


    begin=(-0.4,1.5)
    end=(0.4,1.5)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    begin=(-0.4,-1.5)
    end=(0.4,-1.5)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    begin=(0.4,0.4)
    end=(-0.4,-0.4)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    return None

def boxPlotGrapher(classifiers,borders,flag):

    '''
    this function plots the expression of the classifiers in log10 scale using boxplots
    '''

    # 0. defining a variable for recapitulating the data borders necessary for mapping samples into new space
    borders[flag]={}
    
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

    # 3. plotting and recovering the info for the mapping samples
    boxPlotPosition=0
    listOfClassifiers=sorted(classifiers,key=classifiers.__getitem__,reverse=True) # ranking the classifiers
    for geneID in listOfClassifiers:
        boxPlotPosition=boxPlotPosition+1
        borders[flag][geneID]={}
        
        # 3.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        # 3.2. actual plotting
        if boxplotPlotting == True:
            bp=matplotlib.pyplot.boxplot([logx],positions=[boxPlotPosition],patch_artist=True)
            setBoxColors(bp,'orange')
            bp=matplotlib.pyplot.boxplot([logy],positions=[boxPlotPosition],patch_artist=True)
            setBoxColors(bp,'darkgreen')

        # 3.3. saving the info for the mapping samples
        xa=numpy.min(logx); xb=numpy.median(logx); xc=numpy.max(logx); sdx=numpy.std(logx)
        ya=numpy.min(logy); yb=numpy.median(logy); yc=numpy.max(logy); sdy=numpy.std(logy)
        if xb > yb:
            center=((xa-yc)/2.)+yc
        else:
            center=((ya-xc)/2.)+xc
        # incorporating the data
        #! w=weightRankCalculator(geneID,listOfClassifiers)
        w=weightNMLCalculator(geneID,flag)
        if flag == 'light':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        elif flag == 'growth':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        else:
            print 'error handling flags from boxPlotGrapher (bis)'
            sys.exit()
        
    # 3.3. closing the figure
    if boxplotPlotting == True:
        matplotlib.pyplot.xlim([0,boxPlotPosition+1])
        matplotlib.pyplot.ylim([-0.2,5.])
        theXticks=range(boxPlotPosition)
        theXticksPosition=[element+1 for element in theXticks]
        matplotlib.pyplot.xticks(theXticksPosition,listOfClassifiers,rotation=-90,fontsize=2)
        matplotlib.pyplot.ylabel('log10 FPKM')
        matplotlib.pyplot.tight_layout(pad=2.5)
        #! matplotlib.pyplot.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.33)
        matplotlib.pyplot.savefig('boxplots_%s.pdf'%flag)
        matplotlib.pyplot.clf()

    return borders

def classifiersFilter(classifiers,flag):

    '''
    this function selects the classifiers that have at least separation between max and min values into the two cases in log space
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

def classifiersWriter(selected,flag):

    '''
    this function writes the gene descriptors into a file
    '''

    # f.1. differentiating phases of descriptors
    if flag == 'light':
        f=open('descriptors/epoch0_lightDescriptors_AM_ranked.txt','w')
        g=open('descriptors/epoch0_lightDescriptors_PM_ranked.txt','w')
    elif flag == 'growth':
        f=open('descriptors/epoch0_growthDescriptors_exp_ranked.txt','w')
        g=open('descriptors/epoch0_growthDescriptors_sta_ranked.txt','w')

    sortedList=sorted(selected,key=selected.__getitem__,reverse=True)
    for element in sortedList:
        m=borders[flag][element][1]
        n=borders[flag][element][4]
        if m > n:
            f.write('%s\n'%element)
        else:
            g.write('%s\n'%element)

    f.close()
    g.close()

    # f.2. adding the relative importance of each descriptor
    if flag == 'light':
        f=open('descriptors/epoch0_lightDescriptors_ranked.txt','w')
    elif flag == 'growth':
        f=open('descriptors/epoch0_growthDescriptors_ranked.txt','w')

    accW=0.
    rank=0
    f.write('#rank\tgeneID\tW\taccW\n')
    for element in sortedList:
        rank=rank+1
        localW=borders[flag][element][-3]
        accW=accW+localW
        f.write('%s\t%s\t%.4f\t%.4f\n'%(rank,element,localW,accW))
    f.close()

    return None

def distanceCalculator(sampleID):

    '''
    this function calculates a single value measure of misregulation
    '''

    diurnalLoad=loadCalculator(sampleID,'light')
    growthLoad=loadCalculator(sampleID,'growth')

    if metaData[sampleID]['light'] == 'AM':
        sepx=diurnalLoad
    elif metaData[sampleID]['light'] == 'PM':
        sepx=-diurnalLoad

    if metaData[sampleID]['growth'] == 'exp':
        sepy=growthLoad
    elif metaData[sampleID]['growth'] == 'sta':
        sepy=-growthLoad

    value=0.5*sepx+0.5*sepy
    #! value=numpy.sqrt(sepx**2+sepy**2)
    value=sepx

    #try making each boxplot being one, total 3

    return value

def ellipseSizeCalculator(flag1,flag2):

    '''
    this function calculates size of ellipsoids
    '''

    if flag1 == 'light':
        rankedDimensions=sorted(lightFilteredClassifiers,key=lightFilteredClassifiers.__getitem__,reverse=True)
    elif flag1 == 'growth':
        rankedDimensions=sorted(growthFilteredClassifiers,key=growthFilteredClassifiers.__getitem__,reverse=True)

    averageDispersion=0.
    for transcript in rankedDimensions:

        if flag2 == 'AM' or flag2 == 'exp':
        
            mean=borders[flag1][transcript][1]
            deviation=borders[flag1][transcript][-2]
            top=borders[flag1][transcript][2]
            weight=borders[flag1][transcript][-3]
            
        elif flag2 == 'PM' or flag2 == 'sta':
            
            mean=borders[flag1][transcript][4]
            deviation=borders[flag1][transcript][-1]
            top=borders[flag1][transcript][5]
            weight=borders[flag1][transcript][-3]
     
        
        stretch=top-mean
        if stretch != 0.:
            position=0.5*(deviation*1.96)/stretch
        else:
            position=0.
        value=position*weight
        averageDispersion=averageDispersion+value

    return averageDispersion

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

def loadCalculator(sampleID,flag):

    '''
    this function computes the value of the sample in the new axis 
    '''

    averageLoad=0.
    
    for classifier in borders[flag].keys():
        s=numpy.log10(expression[classifier][sampleID]+1.) # sample expression level

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
                    stretch=abs(NML-f)
                    value=-(NML-s)/stretch
            else:
                if s > NML:
                    stretch=d-NML
                    value=-(s-NML)/stretch
                else:
                    stretch=NML-c
                    value=(NML-s)/stretch
                        
        # dealing with values within previously observed
        else:
            if positive == True: # it is a light sample
                # assuming light boxplot is above
                if b > e:
                    if s > b:
                        stretch=c-b
                        value=1.5+(0.5*(s-b))/stretch
                    else:
                        stretch=b-a
                        value=1.+(0.5*(s-a))/stretch
                # assuming light boxplot is below
                else:
                    if s > b:
                        stretch=c-b
                        value=1.5-(0.5*(s-b))/stretch
                    else:
                        stretch=b-a
                        if stretch != 0.:
                            value=1.5+(0.5*(b-s))/stretch
                        else:
                            value=2.
            else: # it is a dark sample
                # assuming light boxtplot is above
                if b > e:
                    if s > e:
                        stretch=f-e
                        value=-1.5+(0.5*(s-e))/stretch
                    else:
                        stretch=e-d
                        value=-1.5-(0.5*(e-s))/stretch
                        
                # assuming light boxplot is below
                else:
                    if s > e:
                        stretch=f-e
                        value=-1.5-(0.5*(s-e))/stretch
                    else:
                        stretch=e-d
                        value=-1.5+(0.5*(e-s))/stretch

        # weighting the value
        term=w*value
        averageLoad=averageLoad+term

    return averageLoad

def expressionRetriever(classifiers,flag,condition):

    '''
    this function returns the expression values for a set of genes under a specific condition
    '''

    # 1. selecting the acceptable samples
    selectedSamples=[]
    for sampleID in metaData.keys():
        if metaData[sampleID]['epoch'] == 0:
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

                sampleID=vector[6]
                if vector[0] != '':
                    epoch=int(vector[0])
                growth=vector[1]
                light=vector[2]
                co2=int(vector[3].replace(',',''))
                replicate=vector[4]
                preCollapse=bool(int(vector[5]))

                metaData[sampleID]={}
                metaData[sampleID]['growth']=growth
                metaData[sampleID]['epoch']=epoch
                metaData[sampleID]['light']=light
                metaData[sampleID]['co2']=co2
                metaData[sampleID]['replicate']=replicate
                metaData[sampleID]['pre-collapse']=preCollapse

    return metaData

def newCoordinateCalculator(sampleID):

    '''
    this function computes the new coordinates of a sample given the borders calculated from the filtered classifiers
    '''

    x=-loadCalculator(sampleID,'light')
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

        theSize=10
        theAlpha=.85

        # defining the color depending on the light/dark
        if metaData[sampleID]['light'] == 'AM':
            if metaData[sampleID]['replicate'] == 'A':
                theColor='chocolate'
            elif metaData[sampleID]['replicate'] == 'B':
                theColor='orange'
            elif metaData[sampleID]['replicate'] == 'C':
                theColor='orangered'
            else:
                print 'error defining AM replicate at newSpaceMapper. exiting...'
                sys.exit()
        elif metaData[sampleID]['light'] == 'PM':
            if metaData[sampleID]['replicate'] == 'A':
                theColor='olive'
            elif metaData[sampleID]['replicate'] == 'B':
                theColor='green'
            elif metaData[sampleID]['replicate'] == 'C':
                theColor='lightseagreen'
            else:
                print 'error defining PM replicate at newSpaceMapper. exiting...'
                sys.exit()
        else:
            print 'error while defining the color from main'
            sys.exit()

        # defining the marker type depending on epoch
        if metaData[sampleID]['epoch'] == 0:
            theMarker='o'
        elif metaData[sampleID]['epoch'] == 1:
            theMarker='s'
        elif metaData[sampleID]['epoch'] == 2:
            theMarker='^'
        else:
            print 'error while defining the marker from main'
            sys.exit()

        # defining the facecolor and markeredgecolor depending on exp/sta
        if metaData[sampleID]['growth'] == 'exp':
            theMFC='None'; theMEC=theColor
        elif metaData[sampleID]['growth'] == 'sta':
            theMFC=theColor; theMEC='None'
        else:
            print 'error while defining the marker from main'
            sys.exit()
      
        ax.plot(x,y,marker=theMarker,mew=1,color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,zorder=10)

        # marking the pre-collapse samples
        if metaData[sampleID]['pre-collapse'] == True and metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['replicate'] == 'B':
            ax.plot(x,y,marker='*',color='black',ms=3,zorder=10)

    # plotting the ellipses
    a=ellipseSizeCalculator('light','AM')
    b=ellipseSizeCalculator('light','PM')
    c=ellipseSizeCalculator('growth','exp')
    d=ellipseSizeCalculator('growth','sta')

    e=matplotlib.patches.Ellipse(xy=(1.5,1.5),width=b,height=c,edgecolor='darkgreen',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    e=matplotlib.patches.Ellipse(xy=(1.5,-1.5),width=b,height=d,edgecolor='darkgreen',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    
    e=matplotlib.patches.Ellipse(xy=(-1.5,1.5),width=a,height=c,edgecolor='orange',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    e=matplotlib.patches.Ellipse(xy=(-1.5,-1.5),width=a,height=d,edgecolor='orange',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    addArrows()

    # setting ranges
    matplotlib.pyplot.xlim([-2.,2.])
    matplotlib.pyplot.ylim([-2.,2.])
    matplotlib.pyplot.xticks([-1.5,1.5],['light','dark'],fontsize=24)
    matplotlib.pyplot.yticks([-1.5,1.5],['late','early'],fontsize=24)
    matplotlib.pyplot.xlabel('diurnal phase',fontsize=24)
    matplotlib.pyplot.ylabel('growth phase',fontsize=24)
    matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.title(flag+' ppm',fontsize=28)

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

def oneDimensionTimeMapper():

    '''
    this function maps samples into 1D space and plots them as a time series
    '''

    epsilon=0.1
    theSize=15
    theAlpha=.6
    
    # f.1. building a structure holding samples as a time series
    co2levels=[300,1000]
    epochs=[0,1,2]
    growths=['exp','sta']
    diurnals=['AM','PM']

    for co2level in co2levels:
        orderedSamples={}
        time=0
        for epoch in epochs:
            for growth in growths:
                for diurnal in diurnals:
                    time=time+1
                    for sampleID in metaData.keys():
                        if metaData[sampleID]['co2'] == co2level and metaData[sampleID]['epoch'] == epoch and metaData[sampleID]['growth'] == growth and metaData[sampleID]['light'] == diurnal:
                            if time in orderedSamples:
                                orderedSamples[time].append(sampleID)
                            else:
                                orderedSamples[time]=[sampleID]
            
        # f.2. computing distances and building trajectories for each co2 condition
        timePoints=orderedSamples.keys()
        timePoints.sort()
        for timePoint in timePoints:
            samples=orderedSamples[timePoint]

            yPoints=[]
            for sampleID in samples:
                
                yPos=distanceCalculator(sampleID)

                # selecting position and color
                if metaData[sampleID]['replicate'] == 'A':
                    xPos=timePoint-epsilon
                    if co2level == 300:
                        theColor='lightblue'
                    else:
                        theColor='salmon'
                elif metaData[sampleID]['replicate'] == 'B':
                    xPos=timePoint
                    if co2level == 300:
                        theColor='blue'
                    else:
                        theColor='red'
                elif metaData[sampleID]['replicate'] == 'C':
                    xPos=timePoint+epsilon
                    if co2level == 300:
                        theColor='darkblue'
                    else:
                        theColor='darkred'
                else:
                    print 'error choosing the replicates from oneDimensionTimeMapper. exiting...'
                    sys.exit()
            
                # actual plotting the dots
                matplotlib.pyplot.plot(xPos,yPos,marker='.',color=theColor,ms=theSize,alpha=theAlpha,mew=0.,zorder=1)
                yPoints.append(yPos)
            # plotting the bars
            average=numpy.mean(yPoints)
            if co2level == 300:
                theColor='blue'
            else:
                theColor='red'
            matplotlib.pyplot.plot([timePoint-epsilon,timePoint+epsilon],[average,average],'-',color=theColor,lw=2,zorder=2)
             
    # f.3. finishing up the figure
    matplotlib.pyplot.plot([1,12],[1.5,1.5],color='black',ls=':')
    matplotlib.pyplot.plot([1,12],[1.,1.],color='magenta',ls='--')
    matplotlib.pyplot.plot([1,12],[-1.,-1.],color='magenta',ls='--')
    matplotlib.pyplot.plot([4.5,4.5],[-1.5,2],color='black',alpha=0.2,ls='--')
    matplotlib.pyplot.plot([8.5,8.5],[-1.5,2],color='black',alpha=0.2,ls='--')
    
    matplotlib.pyplot.xlabel('epoch',fontsize=24)
    matplotlib.pyplot.ylabel('state',fontsize=24)
    #! matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.11)

    matplotlib.pyplot.xticks((2.5,6.5,10.5),('0','1','2'))
    matplotlib.pyplot.yticks((-1.5,0.,1.5),('inverse','misregulation','expected'))

    matplotlib.pyplot.xlim([0,13])
    matplotlib.pyplot.ylim([-1.5,2])
    
    matplotlib.pyplot.savefig('distance.pdf')
    matplotlib.pyplot.clf()

    return None

def weightNMLCalculator(geneID,flag):

    '''
    this function computes the weight of the classifier based on the empty space between the classifiers
    '''

    if flag == 'light':
        sumSpaces=sum(lightFilteredClassifiers.values())
        value=lightFilteredClassifiers[geneID]
    elif flag == 'growth':
        sumSpaces=sum(growthFilteredClassifiers.values())
        value=growthFilteredClassifiers[geneID]

    weight=value/sumSpaces

    return weight

def weightRankCalculator(geneID,classifiers):

    '''
    this function computes the weight of the classifier based on its rank
    '''

    inverseRank=len(classifiers)-classifiers.index(geneID)
    sumOfRanks=sum(numpy.arange(1.,len(classifiers)+1.))
    weight=float(inverseRank)/sumOfRanks

    return weight

### MAIN

# 0. preliminaries
print 'initializing variables...'
# 0.1. user defined variables and paths
cuffdiffDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'

boxplotPlotting=False
time300=numpy.array([1.375,1.625,3.375,3.708333333,5.291666667,5.708333333,7.333333333,7.75])
time1000=numpy.array([1.375,1.625,3.375,3.708333333,5.291666667,5.708333333,7.333333333,7.75,15.45833333,15.79166667,17.41666667,17.79166667])

# 0.2. reading metadata
metaData=metadataReader()

# 0.3. reading expression
expression=expressionReader()

# 1. recover the DET genes. Rank them on q-value and plot the boxplots. select the best ones.
print 'recovering classifiers...'

lightClassifiers=classifiersRetriever('light_epoch0')
print len(lightClassifiers),'light classifiers detected.'

growthClassifiers=classifiersRetriever('growth_epoch0')
print len(growthClassifiers),'growth classifiers detected.'

# 1.2. filtering the classifiers based on separation
print 'filtering classifiers based on separation...'
lightFilteredClassifiers=classifiersFilter(lightClassifiers,'light')
growthFilteredClassifiers=classifiersFilter(growthClassifiers,'growth')
print len(lightFilteredClassifiers),'filtered light classifiers.'
print len(growthFilteredClassifiers),'filtered growth classifiers.'

# 1.3. plotting a boxplots of the best classifiers
print 'plotting expression graphs for classifiers...'
borders={}
borders=boxPlotGrapher(lightFilteredClassifiers,borders,'light')
borders=boxPlotGrapher(growthFilteredClassifiers,borders,'growth')

classifiersWriter(lightFilteredClassifiers,'light')
classifiersWriter(growthFilteredClassifiers,'growth')

# 2. map samples into a new space of dark/light distributed in x:-2:-1/1:2 and stationary/exponential y:-2:-1/1:2
print 'mapping samples into new space...'
newSpaceMapper('300')
newSpaceMapper('1000')

# 3. map samples into a 1D variable
print 'mapping samples into 1D space...'
oneDimensionTimeMapper()

# 4. final message
print '... analysis completed.'
        
