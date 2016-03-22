import scipy, numpy, sys
from scipy.optimize import leastsq
from scipy.stats import linregress

def characteristicParameterFinder(dataStructure,runningFlags):

    resolution = 1000
    maxGrowthRates = []
    uvValues = []
    growthLags = []
    forceManualFlag=runningFlags[0]
    includeNegativeGrowthFlag=runningFlags[1]

    for epocLabel in dataStructure:

        epoc=epocLabel.split('_')[0]
        uvValue=float(epoc)

        print epocLabel,
                
        localTime=numpy.array(dataStructure[epocLabel][0])
        shiftedTime=localTime-min(localTime)
        localCells=dataStructure[epocLabel][1]
        highResolutionTime=numpy.linspace(min(shiftedTime),max(shiftedTime),resolution)

        # if there are two or less values
        if len(localCells) <= 2:
            maxGrowthRate = None
            growthLag = None
            print 'calculated parameters based on two values', maxGrowthRate,growthLag

        # if the last value is lower than the first, provide with the slope of a regression
        elif localCells[0] > localCells[-1]:
            if includeNegativeGrowthFlag == True:
                if forceManualFlag == True:
                    manualMaxGrowthRate, manualGrowthLag = manualGrowthCalculator(shiftedTime, localCells)
                    maxGrowthRate = manualMaxGrowthRate
                    growthLag = manualGrowthLag
                else:
                    slope, intercept, temp0, temp1, temp2 = scipy.stats.linregress(shiftedTime,localCells)
                    maxGrowthRate=slope
                    growthLag = None
            else:
                maxGrowthRate = None
                growthLag = None
            print 'calculated parameters based on a decreasing time series',maxGrowthRate,growthLag

        # calculate the fitted and the manual max growth. Give the manual if the difference is larger than two fold
        else:
            
            fittedTrajectory = dataFitter(shiftedTime,localCells)
            fittedMaxGrowthRate = fittedTrajectory[0][1]
            manualMaxGrowthRate, manualGrowthLag = manualGrowthCalculator(shiftedTime, localCells)
            ratio = fittedMaxGrowthRate/manualMaxGrowthRate
            
            if ratio > 1.5 or forceManualFlag == True:
                maxGrowthRate = manualMaxGrowthRate
                growthLag = manualGrowthLag
                print 'calculated parameters based on discrete data points',maxGrowthRate,growthLag
            else:
                maxGrowthRate = fittedMaxGrowthRate
                b=(peval(highResolutionTime,fittedTrajectory[0]))
                z=numpy.diff(b,n=1)
                thePosition=numpy.where(z == max(z))
                growthLag=highResolutionTime[thePosition][0]
                print 'calculated parameters based on fitted function', maxGrowthRate, growthLag

        print

        if maxGrowthRate == None or growthLag == None:
            pass
        else:
            uvValues.append(uvValue)
            maxGrowthRates.append(maxGrowthRate)
            growthLags.append(growthLag)

    # working with recovery times
    recoveryTimes=[]
    for value in growthLags:
        if value == 0.:
            recovery=0.
        else:
            recovery=numpy.log(1./value)
        recoveryTimes.append(recovery)
         
    return maxGrowthRates,uvValues,growthLags,recoveryTimes

def dataReader(inputFile):

    dataStructure={} # this dictionary has a label for each epoc. For each epoc you will have a list of two lists for time and number of cells

    with open(inputFile,'r') as f:
        # formating the excell format
        data=f.readlines()
        lines=data[0].split('\r')

        for line in lines:
            vector=line.split('\t')
            if vector[1] != 'Hours':
                
                epoc=vector[0]
                time=float(vector[1])/24.
                cells=float(vector[2])
                
                # filling up the variable
                if epoc not in dataStructure:
                    dataStructure[epoc]=[[time],[cells]]
                else:
                    dataStructure[epoc][0].append(time)
                    dataStructure[epoc][1].append(cells)

    return dataStructure

def dataFitter(time,values):

    p0 = [1e6, 1e5, 1e0, 50e3]
    fitted = scipy.optimize.leastsq(residuals,p0,args=(values,time),maxfev=int(1e9),ftol=1e-20)

    print fitted[0]

    return fitted

def logistic(x, A, B, C, D):
    
    """
    modification of the  logistic equation from grofit, journal of statistical software 2010
    A is the carrying capacity
    B, or mu is the max growth rate
    C, or lambda, length of the lag phase
    D, extra parameter, not for the publication, to account for extra flexibility at t=0
    """
    
    num=A
    expfun=numpy.exp(4*B/A*(C-x)+2)
    deno=1+expfun
    
    return D+num/deno

def manualGrowthCalculator(t,x):

    # checking the sizes
    if len(t) != len(x):
        print len(t),len(x)
        print 'manualGrowthCalculator needs vectors of the same size. Exiting...'
        sys.exit()

    # actual calculation
    mus=[]
    middlePoints = []
    for i in range(len(t)-1):
        
        deltat=t[i+1]-t[i]
        deltax=x[i+1]-x[i]
        mu=deltax/deltat
        mus.append(mu)
        middlePoints.append(t[i] + (deltat)/2.)

    lag=middlePoints[mus.index(max(mus))]

    return max(mus), lag

def peval(x, p):
    """Evaluated value at x with current parameters."""

    A,B,C,D = p

    return logistic(x, A, B, C, D)

def residuals(p, y, x):
    """Deviations of data from fitted 4PL curve"""

    A,B,C,D = p
    err = y-logistic(x, A, B, C, D)

    return err
