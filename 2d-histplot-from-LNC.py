# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 13:10:30 2013

@author: dashamstyr
"""

import LNC_tools as ltools
import os, sys
import numpy as np
import matplotlib.pyplot as plt
import slowhist as h2d
from matplotlib import cm


olddir = os.getcwd()

os.chdir('C:\Users\dashamstyr\Dropbox\PhD - General')

filepath = ltools.get_files('Select LNC file', filetype = ('.h5', '*.h5'))

datatypes = ['BR532','PR532']

header, data_dict = ltools.from_HDF(filepath[0], datatypes)

try:
    BRdat = data_dict[datatypes[0]]
except KeyError:
    sys.exit("No BR532 data available!")
    
try:
    PRdat = data_dict[datatypes[1]]
except KeyError:
    sys.exit("No PR532 data available!")
    
altrange = np.arange(150,10000,10)

BRdat = ltools.alt_resample(BRdat,altrange)
PRdat = ltools.alt_resample(PRdat,altrange)

altrange = BRdat.columns.values

copolvals = np.hstack(BRdat.values).astype('float32')
depolvals = np.hstack(PRdat.values).astype('float32')


numbins = 100
depolmin = 0.0
depolmax = 0.5
copolmin = 0.0
copolmax = 5.0

copolhist=h2d.fullhist(copolvals,numbins,copolmin,copolmax,-9999.,-8888.)
depolhist=h2d.fullhist(depolvals,numbins,depolmin,depolmax,-9999.,-8888.)

altOut = h2d.althist(PRdat.values,altrange,numbins,(depolmin,depolmax))

copolOut=h2d.hist2D(copolhist['fullbins'],depolhist['fullbins'],copolhist['numBins'],depolhist['numBins'])

altcounts=altOut['coverage']
copolcounts = copolOut['coverage']

altcounts[altcounts < 1] = 1
copolcounts[copolcounts < 1] = 1

altlogcounts=np.log10(altcounts)
copollogcounts=np.log10(copolcounts)
                  
try:
    os.chdir('Plots')
except WindowsError:
    os.makedirs('Plots')
    os.chdir('Plots')

startdate = BRdat.index[0].strftime("%Y-%m-%d")
enddate = BRdat.index[-1].strftime("%Y-%m-%d")

starttime = BRdat.index[0].strftime("%H")
endtime = BRdat.index[-1].strftime("%H")

if startdate == enddate:
    if starttime == endtime:
        savetime = startdate+'_'+starttime
    else:
        savetime = startdate+'_'+starttime+'-'+endtime
else:
    savetime = startdate+'-'+enddate

fignum = 0

#fignum+=1
#fig=plt.figure(fignum)
#fig.clf()
#the_axis=fig.add_subplot(111)
#the_axis.plot(depolvals,copolvals,'b+')
#the_axis.set_xlabel('depolvals')
#the_axis.set_ylabel('copolvals')
#the_axis.set_title('raw scatterplot')
#fig.savefig('{0}_{1}-{2}m-copoldepolraw.png'.format(savetime,altrange[0],altrange[-1]))
#fig.canvas.draw()

cmap=cm.bone
cmap.set_over('r')
cmap.set_under('b')

fignum+=1
fig=plt.figure(fignum)
fig.clf()
axis=fig.add_subplot(111)
im=axis.pcolormesh(depolhist['centers'],altrange,altlogcounts.T, cmap = cmap)
cb=plt.colorbar(im,extend='both')
title="2-d histogram"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis.set_xlabel('depolvals')
axis.set_ylabel('Altitude [m]')
axis.set_title(title)
fig.savefig('{0}_{1}-{2}m-althist.png'.format(savetime,altrange[0],altrange[-1]))
fig.canvas.draw()

fignum+=1
fig=plt.figure(fignum)
fig.clf()
axis=fig.add_subplot(111)
im=axis.pcolormesh(depolhist['centers'],copolhist['centers'],copollogcounts, cmap = cmap)
cb=plt.colorbar(im,extend='both')
title="2-d histogram"
colorbar="log10(counts)"
the_label=cb.ax.set_ylabel(colorbar,rotation=270)
axis.set_xlabel('Volume Depolarization Ratio')
axis.set_ylabel('Attenuated Backscatter')
axis.set_title(title)
fig.savefig('{0}_{1}-{2}m-copoldepol.png'.format(savetime,altrange[0],altrange[-1]))
fig.canvas.draw()

#fignum+=1
#fig = plt.figure(fignum)
#fig.clf()
#axis = fig.add_subplot(111)
#hist = axis.hist(depolvals, bins = 100)
#fig.savefig('{0}_{1}-{2}m-1Ddepolhist.png'.format(savetime,altrange[0],altrange[-1]))
#fig.canvas.draw()
plt.show()

os.chdir(olddir)