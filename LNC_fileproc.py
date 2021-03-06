import pandas as pan
import numpy as np
import os,sys
import time, datetime
import re
import csv
import LNC_tools as ltools


#----------------------------------------------------------------------------
#Uses tools created in LNC_tools to open all files in a folder and resample
#them to a regular spacing in altitude/date the concatenates them into one
#pandas dataframe and saves it into an .h5 file format with a header
#Created:     July 05, 2012
#Last Update: Jan 19, 2014
#----------------------------------------------------------------------------


def lnc_reader(filepath):
    #---------------------------------------------------------------------------
    #This program opens an ASCII file generated by the LNC program and outputs a
    #Pandas timeseries dataframe with the data and relevant metadata
    #Note: Input file must be either a 1-D temporal average or a
    #"Single Table with Time in the X Axis vs Altitude in the Y Axis"
    #NOT A "Series of consecutively listed profiles"
    #---------------------------------------------------------------------------

    product = []

    #generate a new filename with _proc attached to the end and open both files
    [fname,fext] = filepath.split('.')
    print 'Pre-Processing '+fname
    fout_name = fname+'_proc'
    fout_name = fout_name+'.'+fext
    fin = open(filepath, 'rb')
    fout = open(fout_name, 'w')

    #copy data table to new *_proc file line by line and replace all spaces
    #between data columns with a comma, empty lines have len(2) & are skipped
    #first line is the product designation, preserved as metadata
    bigspace = re.compile('\s\s\s+')
    for line in fin:
        if not product:
            product.append(line)
        elif len(line) == 2:
            continue  
        else:
            line = bigspace.sub(',',line)+'\n'
            fout.write(line)

    #close both files        
    fin.close()
    fout.close()

    #use csv.reader to read processed file into a list of lists
    temp = []
    for row in csv.reader(open(fout_name,'rb'), delimiter=','):
        temp.append(row)

    #convert to numpy array and transpose list to put datetime entries in first
    #column, which will facilitate conversion to pandas timeseries dataframe
    temparray = np.array(temp).T

    #generate pandas dataframe index by parsing strings into datetime objects
    #note: first entry is the word 'Altitude', last entry is an empty space
    indexdat = []
    for i in temparray[1:-1,0]: indexdat.append(datetime.datetime.strptime(i,'%Y/%m/%d %H:%M:%S')) 
    index = pan.Index(indexdat,name = 'Date Time')
    #generate column headers from altitudes (not including the word 'Altitude'
    coldat = np.array(temparray[0,1:],dtype='float')
    columns = pan.Index(coldat,name = temparray[0,0])
    #data for dataframe consists of remaining rows and columns
    data = temparray[1:-1,1:]

    #check data for flags indicating bad results and substitute with NaN

    flags = ['-1.#INF','1.#INF','-1.#IND','1.#IND']

    clean_data = np.copy(data)

    for f in flags: clean_data[data == f] = np.nan

    #convert data to pandas dataframe
    df = pan.DataFrame(clean_data,index=index,columns=columns,dtype='float')

    return df, product     

def alt_resample(df, altrange):
    #takes a pandas dataframe generated by lnc and resamples on regular
    #intervals in altitude and resets the limits of the set
    #note: limits of altrange must be within original limits of altitude data

    print 'Altitude step resampling in progress ...'
    
    x = df.columns
    
    minalt = df.columns[0]
    maxalt = df.columns[-1]
    
    if minalt > altrange[0]:
        altrange = altrange[altrange >= minalt]
        print "WARNING: Minimum altitude reset to {0}".format(altrange[0])
    
    if maxalt < altrange[-1]:
        altrange = altrange[altrange <= maxalt]
        print "WARNING: Maximum altitude reset to {0}".format(altrange[-1])
    
    numrows = np.size(df.index)
    numcols = np.size(altrange)

    newvalues = np.empty([numrows, numcols])
    n = 0

    for row in df.iterrows():
        f = interp1d(x,row[1].values)
        newvalues[n] = f(altrange)
        n += 1

    dfout = pan.DataFrame(data = newvalues, index = df.index,
                          columns = altrange)
    print '... Done!'
    return dfout

def time_resample(df, timestep, timerange = False, s_mode = 'mean'):
    #resamples a pandas dataframe generated by lnc_reader on a regular timestep
    #and optionally limits it to a preset time range
    #timestep must be in timeseries period format: numF where num=step size and
    #F = offset alias.  Ex: H = hours, M = minutes, S = seconds, L = millieconds
    print 'Time step regularization in progress ...'
    if timerange:
        start_time = timerange[0]
        end_time = timerange[1]

        dfout = df[(df.index>=start_time) & (df.index<=end_time)]

    dfout = dfout.resample(timestep, how = s_mode)

    print '... Done!'
    return dfout    
    

def BR_mask(backscatter, data, delta):
    #this function takes a pandas timeseries dataframe representing a table of
    #backscatter ratios and produces a masking dataframe with values of 0 where
    #ratio is identically 1, and 1 elsewhere then applies it to data
    print 'masking data'
    
    def maskfun(x):
        if x <= 1.0+delta:
            x = np.NaN
        else: x = 1
        
        return x
            
    mask = backscatter.applymap(maskfun)
    masked_data = mask*data
    
#    masked_data.replace(0,np.NaN)
    
    print 'Done!'

    return masked_data

#def save_to_HDF(filename, header, df_dict):
#    
#    #Define class header to create columns to hold header data from .mpl binary file
##    class header(tables.IsDescription):
##        location = tables.StringCol(3)
##        dtype = tables.StringCol(6)        
##        timestamp = tables.Time32Col(1) 
##        timestep = tables.Time32Col(1)
##        numbins = tables.UInt32Col(1) #total number of bins per channel
##        bintime = tables.Float32Col(1)  #bin width in seconds
##        minalt = tables.Float32Col(1) #altitude of lidar station in m AMSL
##        
##       
##    with tables.open_file(filename, mode = 'w', title = 'MPL data file') as h5filename:
##        
##        headertbl = h5filename.create_table('/','Header',header,'Ancillary Data')
##          
##        headerdat = headertbl.row
##                  
##        headerdat['location'] = headerin['location']
##        headerdat['dtype'] = headerin['dtype']      
##        headerdat['timestamp'] = headerin['timestamp']
##        headerdat.append()
##        headertbl.flush()
#        
#    store = pan.HDFStore(filename)
#    store['header'] = header
#    
#    for k,v in df_dict.iteritems():
#        store[k] = v
#        
#    store.close()

#def fromHDF(filename, datatypes):
#    
##    store = pan.HDFStore(filename)
#    
#    header = pan.read_hdf(filename,'header')
#    
#    dataout_dict = {}    
#    for dtype in datatypes:
#        dataout_dict[dtype] = pan.read_hdf(filename,dtype)
#        
##    store.close()
#    
#    return header, dataout_dict

if __name__=='__main__':
    olddir = os.getcwd()
    
    os.chdir('E:\CORALNet\ASCII_Files\Smoke2012\Whistler')
    
    newdir = ltools.set_dir('Select Event Folder')
    
    os.chdir(newdir)
    
    PR_type = 'PR532'
    BR_type = 'BR532'
    MSK_type = 'PR532_msk'
    
    #set altitude range and date step sizes
    minalt = 150  # minimum altitude in meters
    maxalt = 15000  # maximum altitude in meters
    altstep = 30 # altitude step in meters
    altrange = np.arange(minalt,maxalt+altstep,altstep)
    
    start = []#datetime.datetime(2013,04,22,00)
    end = []#datetime.datetime(2013,05,04,17)
    timestep = '600S' #seconds
    
    #set buffer around backscatter ratio of 1 for purposes of creating PR mask
    delta = 0.1
    
    #open event file and parse into file types
    files = os.listdir(newdir)
    BRfiles = []
    PRfiles = []
    procfiles = []
    rawfiles = []
    
    
    
    
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
        elif PR_type in f:
            PRfiles.append(f)
    
    #make sure the files are in a common order of ascending date (assuming they're all
    #from the same station
    BRfiles.sort()
    PRfiles.sort()
    
    #first check to make sure the same number of files in each list
    
    if len(BRfiles) != len(PRfiles):
        sys.exit("Error: Mask files don't match data files")
        
    counter = 0
    header = {}
    
    #double check to make sure the mask files match up with the data files
    for br,pr in zip(BRfiles, PRfiles):
        [br_stat,br_date,br_type] = br.split('_')
        [pr_stat,pr_date,pr_type] = pr.split('_')
        print 'Checking data match for %s'%(br_date)
        if br_date == pr_date and br_stat == pr_stat:
            print 'Check!'
            counter += 1
            if counter == 1:            
                header['location'] = br_stat
        else:
            sys.exit("Error: Data products don't match each other")
    
    #open, altitude resample, and concatenate data and mask files
    
    dataprods = []
    counter = 0
    for br in BRfiles:
        br_temp, tempprod = lnc_reader(br)
        br_realt = alt_resample(br_temp,altrange)
    
        if counter == 0:
            br_event = br_realt
            dataprods.append(BR_type)
        else:
            br_event = pan.concat([br_event,br_realt])
        
        counter +=1
    
    counter = 0
    for pr in PRfiles:
        pr_temp, tempprod = lnc_reader(pr)
        pr_realt = alt_resample(pr_temp,altrange)
    
        if counter == 0:
            pr_event = pr_realt
            dataprods.append(PR_type)
        else:
            pr_event = pan.concat([pr_event,pr_realt])
        
        counter+= 1
            
    
    altrange = br_realt.columns.values    
    #sort by index to make certain data is in order then set date ranges to match
    
    br_event = br_event.sort_index()
    pr_event = pr_event.sort_index()
    
    
    if not start:
        start = br_event.index[0]
    
    if not end:
        end = br_event.index[-1]
    
    br_event = time_resample(br_event,timestep, timerange = [start,end])
    pr_event = time_resample(pr_event,timestep,timerange = [start,end])
    
    msk_event = BR_mask(br_event, pr_event, delta)
    
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
    
    filename = 'LNC_{0}_proc.h5'.format(savetime)
    
    data_dict = {BR_type:br_event,PR_type:pr_event,MSK_type:msk_event}
    
    header = pan.Series(header)
    
    save_to_HDF(filename, header, data_dict)
    
    os.chdir(olddir)