import sys,numpy,matplotlib
import matplotlib.pyplot
import library

def checkFit(a,b,c,d,epocLabel):

    figureFile='results/individualPlots/figure_%s.pdf'%epocLabel
    matplotlib.pyplot.clf()

    # checking for growth first
    if d[0] < d[-1]:
        matplotlib.pyplot.plot(a,b,'-k')
    else:
        print '\t WARNING: omitting model fit because data shows negative growth.'
    
    matplotlib.pyplot.plot(c,d,'ok')

    matplotlib.pyplot.xlim([min(c)-0.6,max(c)+0.6])

    matplotlib.pyplot.ylim([0,1e6])
    matplotlib.pyplot.xlabel('time (days)')
    matplotlib.pyplot.ylabel('number of cells')
    matplotlib.pyplot.savefig(figureFile)

    return None

def dataGrapherSingle(dataStructure,figureLabel):

    resolution=1000

    for epocLabel in dataStructure:

        epoc=epocLabel.split('_')[0]
        localTime=numpy.array(dataStructure[epocLabel][0])
        shiftedTime=localTime-min(localTime)

        localCells=dataStructure[epocLabel][1]
        
        highResolutionTime=numpy.linspace(min(shiftedTime),max(shiftedTime),resolution)

        if len(localCells) > 2: # dealing with sets of at least 2 data points
            print figureLabel+'_'+epocLabel,'\t',

            fittedTrajectory=library.dataFitter(shiftedTime,localCells)
            b=library.peval(highResolutionTime,fittedTrajectory[0])

            checkFit(highResolutionTime,b,shiftedTime,localCells,figureLabel+'_'+epocLabel)

    return None

### MAIN

# 1. data reading
data300=library.dataReader('data/300ppmSet3.txt')
data1000=library.dataReader('data/1000ppmSet3.txt')

# 2. fitting the data to sigmoidal function
print 'fitting data for 300 pppm...'
dataGrapherSingle(data300,'300')

print
print 'fitting data for 1000 ppm...'
dataGrapherSingle(data1000,'1000')

print '... graphs completed.'

