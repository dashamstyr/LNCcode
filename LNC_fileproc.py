import pandas as pan
import numpy as np
import os,sys
import LNC_tools as LNC
import time, datetime

#----------------------------------------------------------------------------
#Uses tools created in LNC_tools to open all files in a folder and resample
#them to a regular spacing in altitude/date the concatenates them into one
#pandas dataframe and plots it using LNC_plot
#July 05, 2012
#----------------------------------------------------------------------------

olddir = os.getcwd()

#os.chdir('K:\CORALNet\Data\ASCII Files')

newdir = LNC.set_dir('Select Event Folder')

os.chdir(newdir)

d_type = 'PR532'
m_type = 'BR532'

files = os.listdir(newdir)
maskfiles = []
datafiles = []
procfiles = []
rawfiles = []

#set altitude range and date step sizes

altrange = np.arange(10,10010,10)#meters
timestep = '120S' #seconds

#set buffer around backscatter ratio of 1 for mask

delta = 0.1

#check to see if each file has been processed before and separate processed
#files into a new list

for f in files:
    if '_proc' in f or '.pickle' in f:
        procfiles.append(f)
    elif '.txt' in f:
        rawfiles.append(f)

#search through list of files to separate fields to be used as a mask from those
#with data to be plotted
#initially, mask files are designated BR1064 for 1064nm Backscatter Ratio

for f in rawfiles:
    if m_type in f:
        maskfiles.append(f)
    elif d_type in f:
        datafiles.append(f)

#make sure the files are in a common order of ascending date (assuming they're all
#from the same station
maskfiles.sort()
datafiles.sort()

#first check to make sure the same number of files in each list

if len(maskfiles) != len(datafiles):
    sys.exit("Error: Mask files don't match data files")
    
counter = 0
dataheader = {}
maskheader = {}
datamaskheader = {}

#double check to make sure the mask files match up with the data files
for d,m in zip(datafiles, maskfiles):
    [d_stat,d_date,d_type] = d.split('_')
    [m_stat,m_date,m_type] = m.split('_')
    print 'Checking mask/data match for %s'%(d_date)
    if d_date == m_date and d_stat == m_stat:
        print 'Check!'
        counter += 1
        if counter == 1:            
            dataheader['location'] = d_stat
            dataheader['dtype'] = d_type
            maskheader['location'] = m_stat
            maskheader['dtype'] = m_type
            datamaskheader['location'] = d_stat
            datamaskheader['dtype'] = dataheader['dtype']+'msk'
        continue
    else:
        sys.exit("Error: Mask files don't match data files")

#open, altitude resample, and concatenate data and mask files

for d,m in zip(datafiles, maskfiles):
    d_temp, data_prod = LNC.lnc_reader(d)
    d_realt = LNC.alt_resample(d_temp,altrange)

    try:
        d_event = pan.concat([d_event,d_realt])
    except NameError:
        d_event = d_realt

    m_temp, data_prod = LNC.lnc_reader(m)
    m_realt = LNC.alt_resample(m_temp,altrange)

    try:
        m_event = pan.concat([m_event,m_realt])
    except NameError:
        m_event = m_realt

    
#sort by index to make certain data is in order then set date ranges to match

d_event = d_event.sort_index()
m_event = m_event.sort_index()

start = m_event.index[0]
end = m_event.index[-1]

d_event = LNC.time_resample(d_event,timestep, timerange = [start,end])
m_event = LNC.time_resample(m_event,timestep,timerange = [start,end])

dt = d_event.index[0].to_pydatetime()
dataheader['timestamp'] = time.mktime(dt.timetuple())
maskheader['timestamp'] = dataheader['timestamp']
datamaskheader['timestamp'] = dataheader['timestamp']

dfmask = LNC.BR_mask(m_event,d_event, delta)

#d_filename = datafiles[0].split('.')[0]+'-'+datafiles[-1].split('.')[0]
#d_event.save(d_filename+'.pickle')
#
#m_filename = maskfiles[0].split('.')[0]+'-'+maskfiles[-1].split('.')[0]
#m_event.save(m_filename+'.pickle')
#
#dfmask.save(d_filename+'_masked.pickle')

d_filename = datafiles[0].split('.')[0]+'-'+datafiles[-1].split('.')[0]+'.h5'
m_filename = maskfiles[0].split('.')[0]+'-'+maskfiles[-1].split('.')[0]+'.h5'
dmask_filename = maskfiles[0].split('.')[0]+'-'+maskfiles[-1].split('.')[0]+'masked.h5'

LNC.save_to_HDF(d_filename, d_event, dataheader)
LNC.save_to_HDF(m_filename, m_event, maskheader)
LNC.save_to_HDF(dmask_filename, dfmask, datamaskheader)
