
def set_dir(titlestring):
    from Tkinter import Tk
    import tkFileDialog
     
    # Make a top-level instance and hide since it is ugly and big.
    root = Tk()
    root.withdraw()
    
    # Make it almost invisible - no decorations, 0 size, top left corner.
    root.overrideredirect(True)
    root.geometry('0x0+0+0')
#    
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    root.deiconify()
    root.attributes("-topmost",1)
    root.focus_force()
    
    file_path = tkFileDialog.askdirectory(parent=root,title=titlestring)
     
    if file_path != "":
       return str(file_path)
     
    else:
       print "you didn't open anything!"
    
    # Get rid of the top-level instance once to make it actually invisible.
    root.destroy() 
     

       
    
     
def get_files(titlestring,filetype = ('.txt','*.txt')):
    from Tkinter import Tk
    import tkFileDialog
    import re
     
     
    # Make a top-level instance and hide since it is ugly and big.
    root = Tk()
    root.withdraw()
    
    # Make it almost invisible - no decorations, 0 size, top left corner.
    root.overrideredirect(True)
    root.geometry('0x0+0+0')
#    
    # Show window again and lift it to top so it can get focus,
    # otherwise dialogs will end up behind the terminal.
    root.deiconify()
    root.attributes("-topmost",1)
    root.focus_force()

    filenames = []
     
    filenames = tkFileDialog.askopenfilename(title=titlestring, filetypes=[filetype,("All files",".")],multiple='True')
    
    #do nothing if already a python list
    if filenames == "": 
        print "You didn't open anything!"  
        return
        
    if isinstance(filenames,list): return filenames

    #http://docs.python.org/library/re.html
    #the re should match: {text and white space in brackets} AND anynonwhitespacetokens
    #*? is a non-greedy match for any character sequence
    #\S is non white space

    #split filenames string up into a proper python list
    result = re.findall("{.*?}|\S+",filenames)

    #remove any {} characters from the start and end of the file names
    result = [ re.sub("^{|}$","",i) for i in result ]     
    return result

    
    root.destroy()

def lnc_reader(filepath):
    #---------------------------------------------------------------------------
    #This program opens an ASCII file generated by the LNC program and outputs a
    #Pandas timeseries dataframe with the data and relevant metadata
    #Note: Input file must be either a 1-D temporal average or a
    #"Single Table with Time in the X Axis vs Altitude in the Y Axis"
    #NOT A "Series of consecutively listed profiles"
    #---------------------------------------------------------------------------
    import pandas as pan
    import numpy as np
    import re
    import csv
    from dateutil.parser import parse

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
    for i in temparray[1:-1,0]: indexdat.append(parse(i)) 
    index = pan.Index(indexdat,name = 'Date Time')
    #generate column headers from altitudes (not including the word 'Altitude'
    coldat = np.array(temparray[0,1:],dtype='float')
    columns = pan.Index(coldat,name = temparray[0,0])
    #data for dataframe consists of remaining rows and columns
    data = temparray[1:-1,1:]

    #check data for flags indicating bad results and substitute with NaN

    flags = ['-1.#INF','1.#INF','-1.#IND','1.#IND']

    clean_data = np.copy(data)

    for f in flags: clean_data[data == f] = 'NaN'

    #convert data to pandas dataframe
    df = pan.DataFrame(clean_data,index=index,columns=columns,dtype='float')

    return df, product     

def alt_resample(df, altrange):
    #takes a pandas dataframe generated by lnc and resamples on regular
    #intervals in altitude and resets the limits of the set
    #note: limits of altrange must be within original limits of altitude data
    
    import numpy as np
    import pandas as pan       
    from scipy.interpolate import interp1d
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
    mask = backscatter.applymap(lambda x: not x <= 1+delta)

    masked_data = mask*data
    
    print 'Done!'

    return masked_data

def save_to_HDF(filename, df, headerin):
    import pandas as pan
    import tables
    import time
    
    #Define class header to create columns to hold header data from .mpl binary file
    class header(tables.IsDescription):
        location = tables.StringCol(3)
        dtype = tables.StringCol(6)        
        timestamp = tables.Time32Col(1) 
        timestep = tables.Time32Col(1)
        numbins = tables.UInt32Col(1) #total number of bins per channel
        bintime = tables.Float32Col(1)  #bin width in seconds
        minalt = tables.Float32Col(1) #altitude of lidar station in m AMSL
        
       
    with tables.open_file(filename, mode = 'w', title = 'MPL data file') as h5filename:
        
        headertbl = h5filename.create_table('/','Header',header,'Ancillary Data')
          
        headerdat = headertbl.row
                  
        headerdat['location'] = headerin['location']
        headerdat['dtype'] = headerin['dtype']      
        headerdat['timestamp'] = headerin['timestamp']
        headerdat.append()
        headertbl.flush()
        
    store = pan.HDFStore(filename)
    store['data'] = df
    store.close()

def from_HDF(filename):
    import pandas as pan
    import tables
    import time, datetime
    
    with tables.openFile(filename,'r+') as h5file: 
        try: 
            table = h5file.root.Header
        except 'tables.exceptions.NoSuchNodeError':
            print 'This file is in the wrong format'
            return
    
        header = {}
        for h in table.iterrows():
            header['location'] = h['location']
            header['dtype'] = h['dtype']      
            header['timestamp'] = datetime.datetime.fromtimestamp(h['timestamp'])
 
        data = pan.read_hdf(filename,'data')

    return data, header

if __name__=='__main__':
    import pandas as pan
    import numpy as np
    import datetime as dt

    delta = 0.1

    BR_filepath = 'C:\Users\dashamstyr\Documents\CORALNet\ASCII Files\UBC_July_2012\UBC_07062012_BR1064.txt'
    data_filepath = 'C:\Users\dashamstyr\Documents\CORALNet\ASCII Files\UBC_July_2012\UBC_07062012_PR532.txt'
    maskout, maskprod = lnc_reader(BR_filepath)
    ar = np.arange(10,18000,100,dtype='float')
    dataout, dataprod = lnc_reader(data_filepath)
    dfmasked = BR_mask(maskout, dataout, delta)

     
