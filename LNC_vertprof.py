# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 15:15:48 2013

@author: dashamstyr
"""

import os, sys
import numpy as np
import datetime as dt
import bisect
import matplotlib.pyplot as plt
import matplotlib.colors as clr
import pandas as pan
import LNC_tools as LNC
import LNC_plot as LNCplt
import datetime as dt

minalt = 300
maxalt = 8000

os.chdir('C:\Users\dashamstyr\Dropbox\PhD - General\My Papers\Smoke Depol July 2012\CORALNet Data')

filepath = LNC.get_files('Select data file', filetype = ('.h5','*.h5'))
            
if filepath[0] == '{':
    filepath = filepath[1:-1]

[path, filename] = os.path.split(filepath[0])

os.chdir(path)

dtypes = ['BR532','PR532_msk']
header, df_dict = LNC.from_HDF(filename, dtypes)

figcount = 0
for k,v in df_dict.iteritems():
    if k in dtypes:    
        df = v
    
        alt = [c for c in df.columns if c >= minalt and c <= maxalt]    
        daterange = df.index
        
        df = df.loc[:,:alt[-1]]
        
        if minalt != 0:
            df.loc[:,:alt[0]] = 'nan'
        
        
        exact_time = []
        approx_time = []
        
        exact_time.append(dt.datetime(2012,7,7,6,0,0))
        exact_time.append(dt.datetime(2012,7,9,12,0,0))
        exact_time.append(dt.datetime(2012,7,12,8,14,0))
        
        #exact_time.append(dt.datetime(2012,7,17,12,0,0))
        #exact_time.append(dt.datetime(2012,7,18,12,0,0))
        #exact_time.append(dt.datetime(2012,7,20,03,0,0))
        
#        exact_time.append(dt.datetime(2012,8,6,12,0,0))
#        exact_time.append(dt.datetime(2012,8,12,6,50,0))
#        exact_time.append(dt.datetime(2012,8,14,10,52,0))
        
        for ts in exact_time:    
            i = bisect.bisect_left(daterange, ts)
            approx_time.append(min(daterange[max(0, i-1): i+2], key=lambda t: abs(ts - t)))
        
        fig = plt.figure(figcount)
        
        numprof = len(approx_time)
        
        for n in range(numprof): 
            print approx_time[n]
            s = df.loc[approx_time[n]]
            #s = s[s>0]
            yvals = s.index
            ymin = yvals[0]
            ymax = yvals[-1]

            ax = fig.add_subplot(1,numprof,n+1)
            
            if k == dtypes[0]:
                im = ax.plot(s,yvals, linewidth = 4)
                plt.ylim([ymin,ymax])
                plt.xlim([1,12])
                plt.xlabel('Backscatter Ratio',fontsize = 21)
            
            if k == dtypes[1]:
                im = ax.scatter(s,yvals)
                plt.ylim([ymin,ymax])
                plt.xlim([0,0.5])
                plt.xlabel('Depolarization Ratio',fontsize = 21)
            
            plt.yticks(fontsize = 21)
            plt.ylabel('Altitude [m]', fontsize = 21)
            
            for line in ax.yaxis.get_ticklines():
                line.set_color('k')
                line.set_markersize(6)
                line.set_markeredgewidth(2)
                
            for line in ax.xaxis.get_ticklines():
                line.set_color('k')
                line.set_markersize(6)
                line.set_markeredgewidth(2)    
                
            plt.xticks(fontsize = 21)
            plt.title(approx_time[n], fontsize = 21)
            fig.subplots_adjust(wspace = 0.5)
            figcount +=1

plt.show()
    