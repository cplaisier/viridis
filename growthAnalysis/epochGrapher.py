import sys,numpy,matplotlib
import matplotlib.pyplot, scipy.stats
import library

def colorDefiner(epoch):

    if epoch == '0':
        theColor='blue'
    elif epoch == '0.5':
        theColor='red'
    elif epoch == '1':
        theColor='green'
    elif epoch == '1.5':
        theColor='orange'
    else:
        print 'error from colorDefiner. exiting...'
        sys.exit()

    return theColor

def dataGrapherEpochs(dataStructure,figureLabel):

    resolution=1000
    figureFile='results/figure_%s.pdf'%figureLabel

    for epochLabel in dataStructure:


        epoch=epochLabel.split('_')[0]
        
        localTime=numpy.array(dataStructure[epochLabel][0])
        shiftedTime=localTime-min(localTime)
        localCells=dataStructure[epochLabel][1]
        highResolutionTime=numpy.linspace(min(shiftedTime),max(shiftedTime),resolution)

        epochColor=colorDefiner(epoch)

        # plotting the data
        if len(localCells) > 1:
            matplotlib.pyplot.plot(localTime,localCells,'o',color=epochColor,markeredgecolor='None',ms=4)

        # plotting the model if there is growth, otherwise plot a best model straight line
        if len(localCells) <= 2:
            matplotlib.pyplot.plot([localTime[0],localTime[-1]],[localCells[0],localCells[-1]],'-',color=epochColor)
        elif localCells[0] > localCells[-1]:
            slope, intercept, temp0, temp1, temp2 = scipy.stats.linregress(shiftedTime,localCells)
            matplotlib.pyplot.plot([localTime[0],localTime[-1]],[intercept,slope*shiftedTime[-1]+intercept],'-',color=epochColor)
        else:
            fittedTrajectory=library.dataFitter(shiftedTime,localCells)
            b=library.peval(highResolutionTime,fittedTrajectory[0])
            matplotlib.pyplot.plot(highResolutionTime+min(localTime),b,'-',color=epochColor)
    
    matplotlib.pyplot.xlim([-0.5,20])
    matplotlib.pyplot.ylim([-0.5e5,18e5])
    matplotlib.pyplot.xlabel('time (days)')
    matplotlib.pyplot.ylabel('number of cells (x 1e5)')
    matplotlib.pyplot.title('%s ppm'%figureLabel)
    matplotlib.pyplot.yticks((0,2e5,4e5,6e5,8e5,10e5,12e5,14e5,16e5,18e5),('0','2','4','6','8','10','12','14','16','18'))

    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()


    return None

### MAIN

# 1. data reading
data300=library.dataReader('data/300ppmSetsLight.v2.txt')
data1000=library.dataReader('data/1000ppmSetsLight.v2.txt')

# 2. fitting the data to sigmoidal function
print 'fitting data for 300 pppm...'
dataGrapherEpochs(data300,'300')

print
print 'fitting data for 1000 ppm...'
dataGrapherEpochs(data1000,'1000')

print '... graphs completed.'
