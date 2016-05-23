### this script generates two cartoon figures about the state space calculation

import sys,matplotlib,numpy
from matplotlib import pyplot

def setBoxColors(bp,theColor):

    '''
    this function access the elements of a boxplot and colors them appropriately
    '''

    matplotlib.pyplot.setp(bp['boxes'],color=theColor)
    matplotlib.pyplot.setp(bp['caps'],color='None')
    matplotlib.pyplot.setp(bp['whiskers'],color=theColor,ls='-')
    matplotlib.pyplot.setp(bp['fliers'],markeredgecolor=theColor,marker='+')
    matplotlib.pyplot.setp(bp['medians'],color=theColor)    

    return None

# 1. building the figure about samples position

x=[-1.5,1.5,0.,0.,1.5]
y=[1.5,0,0.,1.5,-1.5]
theMarkers=['o','v','D','^','s']

theSize=10
theAlpha=.85
theColor='black'
theMFC=theColor; theMEC='None'

for i in range(len(x)):
    matplotlib.pyplot.plot(x[i],y[i],marker=theMarkers[i],color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=0.)


# setting ranges
matplotlib.pyplot.xlim([-2.,2.])
matplotlib.pyplot.ylim([-2.,2.])
matplotlib.pyplot.xticks([-1.5,1.5],['light','dark'],fontsize=24)
matplotlib.pyplot.yticks([-1.5,1.5],['late','early'],fontsize=24)
matplotlib.pyplot.xlabel('diurnal cycle',fontsize=24)
matplotlib.pyplot.ylabel('growth phase',fontsize=24)

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
matplotlib.pyplot.plot([-1,-1],[-1,1],color='black',alpha=0.7,lw=2.)
matplotlib.pyplot.plot([-1,1],[1,1],color='black',alpha=0.7,lw=2.)
matplotlib.pyplot.plot([1,1],[1,-1],color='black',alpha=0.7,lw=2.)
matplotlib.pyplot.plot([-1,1],[-1,-1],color='black',alpha=0.7,lw=2.)

# aspect
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.axes().set_aspect('equal')

matplotlib.pyplot.savefig('figuresCartoon/cartoonStates.pdf')
matplotlib.pyplot.clf()

    
matplotlib.pyplot.savefig('spaceCartoon.pdf')
matplotlib.pyplot.clf()

# 2. building the figure about the box plots
names=['A','B','C','D']
flags=['diurnal','diurnal','growth','growth']
lowx=numpy.log10([1.,5.,1.,5.])
highx=numpy.log10([10.,50.,10.,50.])
lowy=numpy.log10([100.,500.,100.,500.])
highy=numpy.log10([1000.,5000.,1000.,5000.])
resolution=1e3
pos=[1,2,3,4]
theSize=7
theAlpha=.85
theColor='black'
theMFC=theColor; theMEC='None'

matplotlib.pyplot.figure()

# 1. early light
matplotlib.pyplot.subplot(2,2,1)
boxPlotPosition=0.
for i in range(len(names)):
    boxPlotPosition=boxPlotPosition+1

    x=numpy.random.uniform(lowx[i],highx[i],resolution)
    y=numpy.random.uniform(lowy[i],highy[i],resolution)

    bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'orange')
    else:
        setBoxColors(bp,'#0571b0')
    bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'darkgreen')
    else:
        setBoxColors(bp,'#ca0020')

a=0.5
b=numpy.log10(5.)*0.5+numpy.log10(50.)*0.5
sig=([a,b,a,b])
matplotlib.pyplot.plot(pos,sig,marker='o',color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=1.,ls=':',zorder=10)

matplotlib.pyplot.xlim([0.5,4.5])
matplotlib.pyplot.ylim([-0.2,4.])
theXticks=range(4)
theXticksPosition=[element+1 for element in theXticks]
matplotlib.pyplot.xticks(theXticksPosition,names,fontsize=14)
matplotlib.pyplot.yticks([1,2,3,4])
matplotlib.pyplot.ylabel('log10 FPKM')
matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off')

# 4. late dark
matplotlib.pyplot.subplot(222)
boxPlotPosition=0.
for i in range(len(names)):
    boxPlotPosition=boxPlotPosition+1

    x=numpy.random.uniform(lowx[i],highx[i],resolution)
    y=numpy.random.uniform(lowy[i],highy[i],resolution)

    bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'orange')
    else:
        setBoxColors(bp,'#0571b0')
    bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'darkgreen')
    else:
        setBoxColors(bp,'#ca0020')

a=2.5
b=numpy.log10(500.)*0.5+numpy.log10(5000.)*0.5
sig=([a,b,a,b])
matplotlib.pyplot.plot(pos,sig,marker='s',color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=1.,ls=':',zorder=10)

matplotlib.pyplot.xlim([0.5,4.5])
matplotlib.pyplot.ylim([-0.2,4.])
theXticks=range(4)
theXticksPosition=[element+1 for element in theXticks]
matplotlib.pyplot.xticks(theXticksPosition,names,fontsize=14)
matplotlib.pyplot.yticks([1,2,3,4])
matplotlib.pyplot.ylabel('log10 FPKM')
matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off')

# 3. total mis regulation
matplotlib.pyplot.subplot(223)
boxPlotPosition=0.
for i in range(len(names)):
    boxPlotPosition=boxPlotPosition+1

    x=numpy.random.uniform(lowx[i],highx[i],resolution)
    y=numpy.random.uniform(lowy[i],highy[i],resolution)

    bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'orange')
    else:
        setBoxColors(bp,'#0571b0')
    bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'darkgreen')
    else:
        setBoxColors(bp,'#ca0020')

c=1.5
d=0.5+numpy.log10(50.)
sig=([c,d,c,d])
matplotlib.pyplot.plot(pos,sig,marker='D',color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=1.,ls=':',zorder=10)

matplotlib.pyplot.xlim([0.5,4.5])
matplotlib.pyplot.ylim([-0.2,4.])
theXticks=range(4)
theXticksPosition=[element+1 for element in theXticks]
matplotlib.pyplot.xticks(theXticksPosition,names,fontsize=14)
matplotlib.pyplot.yticks([1,2,3,4])
matplotlib.pyplot.ylabel('log10 FPKM')
matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off')

# 2. partial lack of regulation
matplotlib.pyplot.subplot(224)
boxPlotPosition=0.
for i in range(len(names)):
    boxPlotPosition=boxPlotPosition+1

    x=numpy.random.uniform(lowx[i],highx[i],resolution)
    y=numpy.random.uniform(lowy[i],highy[i],resolution)

    bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'orange')
    else:
        setBoxColors(bp,'#0571b0')
    bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
    if flags[i] == 'diurnal':
        setBoxColors(bp,'darkgreen')
    else:
        setBoxColors(bp,'#ca0020')

a=0.5
b=numpy.log10(5.)*0.5+numpy.log10(50.)*0.5
c=1.5
d=0.5+numpy.log10(50.)
sig1=[c,d,a,b]

a=2.5
b=numpy.log10(500.)*0.5+numpy.log10(5000.)*0.5
sig2=[a,b,c,d]
matplotlib.pyplot.plot(pos,sig1,marker='^',color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=1.,ls=':',zorder=10)
matplotlib.pyplot.plot(pos,sig2,marker='v',color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,lw=1.,ls=':',zorder=10)

matplotlib.pyplot.xlim([0.5,4.5])
matplotlib.pyplot.ylim([-0.2,4.])
theXticks=range(4)
theXticksPosition=[element+1 for element in theXticks]
matplotlib.pyplot.xticks(theXticksPosition,names,fontsize=14)
matplotlib.pyplot.yticks([1,2,3,4])
matplotlib.pyplot.ylabel('log10 FPKM')
matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
matplotlib.pyplot.tick_params(axis='y',which='both',right='off')

# setting the legend
matplotlib.pyplot.errorbar([-1, -1, -1], [2, 3, 1],yerr=0.4, fmt="s",label='light',color='orange',lw=2,mec='orange')
matplotlib.pyplot.errorbar([-1, -1, -1], [2, 3, 1],yerr=0.4, fmt="s",label='dark',color='darkgreen',lw=2,mec='darkgreen')
matplotlib.pyplot.errorbar([-1, -1, -1], [2, 3, 1],yerr=0.4, fmt="s",label='early',color='#0571b0',lw=2,mec='#0571b0')
matplotlib.pyplot.errorbar([-1, -1, -1], [2, 3, 1],yerr=0.4, fmt="s",label='late',color='#ca0020',lw=2,mec='#ca0020')
matplotlib.pyplot.legend(numpoints=1,loc='upper center', bbox_to_anchor=(-0.175, -0.12), ncol=4)
matplotlib.pyplot.tight_layout(pad=3.)

matplotlib.pyplot.savefig('figuresCartoon/cartoonBoxplot.pdf')
matplotlib.pyplot.clf()
