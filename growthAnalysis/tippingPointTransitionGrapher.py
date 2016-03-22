### this script plots the tipping point transition based on max growth and growth lag

import matplotlib,sys,getopt,numpy
import matplotlib.pyplot
import library

def markerDefiner(m):

    if m == 0.:
        theMarker = 'o'
    elif m == 0.5:
        theMarker = 's'
    elif m == 1.:
        theMarker = '^'
    elif m == 1.5:
        theMarker = 'v'
    else:
        print 'error trying to assign colors. exiting...'
        sys.exit()

    return theMarker

### MAIN

# 0. reading input options
forceManualFlag=False
includeNegativeGrowthFlag=False
try:
    opts, args = getopt.getopt(sys.argv[1:], 'fn', ['force-manual', 'negative'])
    for o, a in opts:
        if o in ("-f", "--force-manual"):
            forceManualFlag = True
        elif o in ('-n', '--negative'):
            includeNegativeGrowthFlag = True
except:
    print 'ERROR: only flags admitted are -f [--force-manual] and -n [--negative].'
    sys.exit()
    
runningFlags=[forceManualFlag,includeNegativeGrowthFlag]

# 1. data reading
data300 = library.dataReader('data/300ppmSetsLight.v2.txt')
data1000 = library.dataReader('data/1000ppmSetsLight.v2.txt')

# 2. calculating the max growth rates
print 'fitting data for 300 pppm...'
maxGrowthRates300, uvValues300, growthLag300, recovery300 = library.characteristicParameterFinder(data300,runningFlags)

print
print 'fitting data for 1,000 pppm...'
maxGrowthRates1000, uvValues1000, growthLag1000, recovery1000 = library.characteristicParameterFinder(data1000,runningFlags)


# 3. plotting
print
print 'plotting the figure...'
figureFile='results/figureTPT'
if runningFlags[0] == True:
    figureFile=figureFile+'_forcedManual'
if runningFlags[1] == True:
    figureFile=figureFile+'_withNegativeGrowth'
figureFile=figureFile+'.pdf'

# 4.1 plotting for 300 ppms
print 'plotting for 300 ppms...'

dataStructurex={}
dataStructurey={}

for i in range(len(maxGrowthRates300)):
    x=growthLag300[i]
    y=maxGrowthRates300[i]
    m=uvValues300[i]

    if x != None:
        theMarker = markerDefiner(m)
        matplotlib.pyplot.plot(x, y, theMarker, color='blue', mec='None', mfc='blue', ms=4, mew=1, alpha=0.5)
        
        if theMarker not in dataStructurex.keys():
            dataStructurex[theMarker]=[x]
            dataStructurey[theMarker]=[y]
        else:
            dataStructurex[theMarker].append(x)
            dataStructurey[theMarker].append(y)

average_x=[]; average_y=[]; average_errx=[]; average_erry=[]
for theMarker in dataStructurex:
    average_x.append(numpy.mean(dataStructurex[theMarker]))
    average_y.append(numpy.mean(dataStructurey[theMarker]))
    average_errx.append(numpy.std(dataStructurex[theMarker]))
    average_erry.append(numpy.std(dataStructurey[theMarker]))
    matplotlib.pyplot.errorbar(average_x,average_y,xerr=average_errx,yerr=average_erry,fmt=theMarker,color='blue',mec='None')        

# 4.2 plotting for 1,000 ppms
print 'plotting for 1,000 ppms...'

dataStructurex={}
dataStructurey={}

for i in range(len(maxGrowthRates1000)):
    x=growthLag1000[i]
    y=maxGrowthRates1000[i]
    m=uvValues1000[i]

    if x != None:
        theMarker = markerDefiner(m)
        matplotlib.pyplot.plot(x, y, theMarker, color='red', mec='None', mfc='red', ms=4, mew=1, alpha=0.5)

        if theMarker not in dataStructurex.keys():
            dataStructurex[theMarker]=[x]
            dataStructurey[theMarker]=[y]
        else:
            dataStructurex[theMarker].append(x)
            dataStructurey[theMarker].append(y)

average_x=[]; average_y=[]; average_errx=[]; average_erry=[]
for theMarker in dataStructurex:
    average_x.append(numpy.mean(dataStructurex[theMarker]))
    average_y.append(numpy.mean(dataStructurey[theMarker]))
    average_errx.append(numpy.std(dataStructurex[theMarker]))
    average_erry.append(numpy.std(dataStructurey[theMarker]))
    matplotlib.pyplot.errorbar(average_x,average_y,xerr=average_errx,yerr=average_erry,fmt=theMarker,color='red',mec='None')

matplotlib.pyplot.xlim(0,9)
matplotlib.pyplot.xlabel('Max Growth Time Lag (day)')


matplotlib.pyplot.ylabel('Max Growth (cells x 1e5/day)')
matplotlib.pyplot.yticks((0, 2e5, 4e5, 6e5, 8e5, 10e5, 12e5),('0','2','4','6','8','10','12'))

lowesty=min(maxGrowthRates300+maxGrowthRates1000)
highesty=max(maxGrowthRates300+maxGrowthRates1000)
matplotlib.pyplot.ylim(lowesty-0.1*highesty,highesty+0.1*highesty)

matplotlib.pyplot.plot(-1,-1,'-',lw=4, color='blue',label='300 ppm')
matplotlib.pyplot.plot(-1,-1,'-',lw=4,color='red',label='1,000 ppm')
matplotlib.pyplot.legend(loc=1,frameon=False)
matplotlib.pyplot.text(3,7e5,'0.0 r.u. UV',color='black')
matplotlib.pyplot.text(4,3.5e5,'0.5 r.u. UV',color='black')
matplotlib.pyplot.text(6.8,3e5,'1.0 r.u. UV',color='black')

matplotlib.pyplot.savefig(figureFile)

print '... graphs completed.'
