import pandas as pan
import numpy as np
import os,sys
import LNC_tools as LNC
import time, datetime

#----------------------------------------------------------------------------
#Uses tools created in LNC_tools to open all files in a folder and resample
#them to a regular spacing in altitude/date the concatenates them into one
#pandas dataframe and plots it using LNC_plot
#This version only handles backscatter ratio, not polarization
#October 07, 2013
#----------------------------------------------------------------------------

olddir = os.getcwd()

#os.chdir('K:\CORALNet\Data\ASCII Files')

newdir = LNC.set_dir('Select Event Folder')

os.chdir(newdir)

BR_type = 'BR1064'

files = os.listdir(newdir)
BRfiles = []
PRfiles = []
procfiles = []
rawfiles = []

#set altitude range and date step sizes

altrange = np.arange(10,5010,10)#meters
timestep = '120S' #seconds

#set buffer around backscatter ratio of 1 for mask

#delta = 0.1

#check to see if each file has been processed before and separate processed
#files into a new list

for f in files:
    if '_proc' in f or '.h5' in f:
        procfiles.append(f)
    elif '.txt' in f:
        rawfiles.append(f)

#search through list of files to separate fields to be used as a mask from those
#with data to be plotted
#initially, mask files are designated BR1064 for 1064nm Backscatter Ratio

for f in rawfiles:
    if BR_type in f:
        BRfiles.append(f)

#make sure the files are in a common order of ascending date (assuming they're all
#from the same station
BRfiles.sort()
    
counter = 0
header = {}

#open, altitude resample, and concatenate data and mask files

dataprods = []
counter = 0
for br in BRfiles:
    br_temp, tempprod = LNC.lnc_reader(br)
    br_realt = LNC.alt_resample(br_temp,altrange)

    if counter == 0:
        br_event = br_realt
        dataprods.append(BR_type)
    else:
        br_event = pan.concat([br_event,br_realt])
    
    counter +=1
        

altrange = br_realt.columns.values    
#sort by index to make certain data is in order then set date ranges to match

br_event = br_event.sort_index()

start = br_event.index[0]
end = br_event.index[-1]

br_event = LNC.time_resample(br_event,timestep, timerange = [start,end])

dt = br_event.index[0].to_pydatetime()
header['timestamp'] = dt
header['timestep'] = timestep
header['dataprods'] = dataprods
header['numbins'] = len(altrange)
header['numprofs'] = len(br_event.index)
header['minalt'] = altrange[0]
header['maxalt'] = altrange[-1]


startdate = start.strftime("%Y-%m-%d")
enddate = end.strftime("%Y-%m-%d")

starttime = start.strftime("%H")
endtime = end.strftime("%H")

if startdate == enddate:
    if starttime == endtime:
        savetime = startdate+'_'+starttime
    else:
        savetime = startdate+'_'+starttime+'-'+endtime
else:
    savetime = startdate+'-'+enddate

filename = 'LNC_{0}.h5'.format(savetime)

data_dict = {BR_type:br_event}

header = pan.Series(header)

LNC.save_to_HDF(filename, header, data_dict)
#LNC.save_to_HDF(m_filename, m_event, maskheader)
#LNC.save_to_HDF(dmask_filename, dfmask, datamaskheader)
