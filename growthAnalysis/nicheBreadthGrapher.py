### this script plots the niche breadth increase due to max growth and UV

import sys, numpy, scipy, matplotlib, getopt
import matplotlib.pyplot
import library

### MAIN

# 0. reading input options
forceManualFlag=False
includeNegativeGrowthFlag=False
try:
    opts, args = getopt.getopt(sys.argv[1:], 'fn', ['force-manual', 'negative'])
    for o, a in opts:
        if o in ('-f', '--force-manual'):
            forceManualFlag = True
        elif o in ('-n', '--negative'):
            includeNegativeGrowthFlag = True
except:
    print 'ERROR: only flags admitted are -f [--force-manual] and -n [--negative].'
    sys.exit()
    
runningFlags=[forceManualFlag,includeNegativeGrowthFlag]

# 1. data reading
data300=library.dataReader('data/300ppmSetsLight.v2.txt')
data1000=library.dataReader('data/1000ppmSetsLight.v2.txt')

# 2. calculating the max growth rates
print 'fitting data for 300 pppm...'
maxGrowthRates300, uvValues300, growthLag300, recovery300 = library.characteristicParameterFinder(data300,runningFlags)

print
print 'fitting data for 1,000 pppm...'
maxGrowthRates1000, uvValues1000, growthLag1000, recovery1000 = library.characteristicParameterFinder(data1000,runningFlags)

# 3. plotting
print
print 'plotting the figure...'
figureFile='results/figureNB'
if runningFlags[0] == True:
    figureFile=figureFile+'_forcedManual'
if runningFlags[1] == True:
    figureFile=figureFile+'_withNegativeGrowth'
figureFile=figureFile+'.pdf'

shift=0.01
# 300
matplotlib.pyplot.plot(numpy.array(uvValues300)-shift, maxGrowthRates300, '.', color='blue', mec='None', mfc='blue', ms=8, mew=1, alpha=0.5)

dataStructure={}
for i in range(len(uvValues300)):
    if uvValues300[i] not in dataStructure.keys():
        dataStructure[uvValues300[i]]=[maxGrowthRates300[i]]
    else:
        dataStructure[uvValues300[i]].append(maxGrowthRates300[i])

average_x=[]; average_y=[]; average_err=[]
for element in dataStructure:
    average_x.append(element)
    average_y.append(numpy.mean(dataStructure[element]))
    average_err.append(numpy.std(dataStructure[element]))

matplotlib.pyplot.errorbar(numpy.array(average_x)+shift, average_y, yerr=average_err, fmt='o', color='blue', mec='None')

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(uvValues300, maxGrowthRates300)
y300 = slope*numpy.array(uvValues300) + intercept
matplotlib.pyplot.plot(uvValues300, y300, color='blue', lw=2, label='300 ppm')
matplotlib.pyplot.text(0.1,2e5,'r=%.3f, p=%.3f'%(r_value,p_value),color='blue')

# 1,000
matplotlib.pyplot.plot(numpy.array(uvValues1000)-shift, maxGrowthRates1000, '.', color='red', mec='None', mfc='red', ms=8, mew=1, alpha=0.5)

dataStructure={}
for i in range(len(uvValues1000)):
    if uvValues1000[i] not in dataStructure.keys():
        dataStructure[uvValues1000[i]]=[maxGrowthRates1000[i]]
    else:
        dataStructure[uvValues1000[i]].append(maxGrowthRates1000[i])

average_x=[]; average_y=[]; average_err=[]
for element in dataStructure:
    average_x.append(element)
    average_y.append(numpy.mean(dataStructure[element]))
    average_err.append(numpy.std(dataStructure[element]))

matplotlib.pyplot.errorbar(numpy.array(average_x)+shift, average_y, yerr=average_err, fmt='o', color='red', mec='None')

slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(uvValues1000, maxGrowthRates1000)
y1000 = slope*numpy.array(uvValues1000) + intercept
matplotlib.pyplot.plot(uvValues1000, y1000, color='red', lw=2, label='1,000 ppm')
matplotlib.pyplot.text(0.1,8e5,'r=%.3f, p=%.3f'%(r_value,p_value),color='red')

low_xlim=min(uvValues1000)-0.1*max(uvValues1000)
up_xlim=max(uvValues1000)+0.1*max(uvValues1000)
matplotlib.pyplot.xlim(low_xlim,up_xlim)

minValue=min(maxGrowthRates300+maxGrowthRates1000+list(y300))
maxValue=max(maxGrowthRates300+maxGrowthRates1000)
low_ylim=minValue-0.1*maxValue
up_ylim=maxValue+0.1*maxValue
matplotlib.pyplot.ylim(low_ylim,up_ylim)

matplotlib.pyplot.xlabel('UV (ru)')
matplotlib.pyplot.ylabel('max growth (cells x 1e5/day)')
matplotlib.pyplot.yticks((0, 2e5, 4e5, 6e5, 8e5, 10e5, 12e5),('0','2','4','6','8','10','12'))
matplotlib.pyplot.legend(loc=1,frameon=False)
matplotlib.pyplot.savefig(figureFile)

print '... graphs completed.'
