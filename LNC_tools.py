import os,sys,site
home=os.environ['homepath']

from Tkinter import Tk
import tkFileDialog
import re
import numpy as np
import array, struct
import pandas as pan
import csv
from dateutil.parser import parse
import datetime
from scipy import constants as const    
from copy import deepcopy
from scipy.interpolate import interp1d
from scipy.ndimage.filters import generic_filter as genfilt
import matplotlib.pyplot as plt
from collections import OrderedDict 
import h5py 


class LNC:
    """
    This is a class type generated by unpacking a binary file generated by
    the CORALNet lidar
    
    it includes a 
    
    """
        
    def __init__(self,filename=[]):
        
        self.BR = None #slot for backscatter ratio dataframe
        self.PR = None  #slot for rdepolarization ratio dataframe
        self.MSK = None   #slot for data mask
        self.sigma = None   #slot for standard deviation data
        self.SNR = None     #slot for SNR data
        self.backscatter = None #slot for corrected backscatter array
        self.extinction = None  #slot for extinction array
        self.scenepanel = None  #slot for panel containing scene analysis features
        
    def frompickle(self,BRfile,PRfile,MSKfile,altcor=True):
        self.BR=pan.read_pickle(BRfile)
        self.PR=pan.read_pickle(PRfile)
        self.MSK=pan.read_pickle(MSKfile)
        
        if altcor:
            #if this is true, change altitude units from meters to kilometers       
            for df in [self.BR,self.PR,self.MSK]:
                df.rename(columns=lambda x: x/1000.0,inplace=True)
        return self
        
    def fromHDF(self, filename,verbose = False):
                
        BRdat = pan.read_hdf(filename,'BR')
        PRdat = pan.read_hdf(filename,'PR')
        MSKdat = pan.read_hdf(filename,'MSK')
            
        self.BR = BRdat
        self.PR = PRdat
        self.MSK = MSKdat
                
        
        try:
            backdat=pan.read_hdf(filename,'backscatter')
            self.backscatter=backdat
        except KeyError:
            if verbose:
                print "Warning: No Backscatter file"
        
        try:
            extdat=pan.read_hdf(filename,'extinction')
            self.extinction=extdat
        except KeyError:
            if verbose:
                print "Warning: No Extinction file"
        
        try:
            scenedat=pan.read_hdf(filename,'scenepanel')
            self.scenepanel=scenedat
        except KeyError:
            if verbose:
                print "Warning: No Scene Analysis file"
                
        sigmadict={}
        try:            
            tempsigma_BR=pan.read_hdf(filename,'sigma_BR')
            sigmadict['BR']=tempsigma_BR
        except KeyError:
            if verbose:
                print "Warning: No sigma-BR file"
        try:            
            tempsigma_PR=pan.read_hdf(filename,'sigma_PR')
            sigmadict['PR']=tempsigma_PR
        except KeyError:
            if verbose:
                print "Warning: No sigma-PR  file"

        try:            
            tempsigma_back=pan.read_hdf(filename,'sigma_backscatter')
            sigmadict['backscatter']=tempsigma_back
        except KeyError:
            if verbose:
                print "Warning: No sigma-backscatter file"

        try:            
            tempsigma_ext=pan.read_hdf(filename,'sigma_extinction')
            sigmadict['extinction']=tempsigma_ext
        except KeyError:
            if verbose:
                print "Warning: No sigma-extinction file"
        
        if sigmadict:
            self.sigma=sigmadict
            
        SNRdict={}
        try:            
            tempSNR_BR=pan.read_hdf(filename,'SNR_BR')
            SNRdict['BR']=tempSNR_BR
        except KeyError:
            if verbose:
                print "Warning: No SNR-BR file"
        try:            
            tempSNR_PR=pan.read_hdf(filename,'SNR_PR')
            SNRdict['PR']=tempSNR_PR
        except KeyError:
            if verbose:
                print "Warning: No SNR-PR  file"

        try:            
            tempSNR_back=pan.read_hdf(filename,'SNR_backscatter')
            SNRdict['backscatter']=tempSNR_back
        except KeyError:
            if verbose:
                print "Warning: No SNR-backscatter file"

        try:            
            tempSNR_ext=pan.read_hdf(filename,'SNR_extinction')
            SNRdict['extinction']=tempSNR_ext
        except KeyError:
            if verbose:
                print "Warning: No SNR-extinction file"
                
        if SNRdict:
            self.SNR=SNRdict
            
        return self
        
    def save_to_HDF(self, filename):
        
        store = pan.HDFStore(filename)
        
        df_BR = self.BR
        store['BR'] = df_BR
        df_PR = self.PR
        store['PR'] = df_PR
        df_MSK = self.MSK
        store['MSK'] = df_MSK
        
        if self.backscatter is not None:
            df_backscatter = self.backscatter
            store['backscatter'] = df_backscatter
        
        if self.extinction is not None:
            df_extinction = self.extinction
            store['extinction'] = df_extinction
        
        if self.scenepanel is not None:
            scenepanel = self.scenepanel
            store['scenepanel'] = scenepanel
            
        if self.SNR is not None:
            for k,v in self.SNR.iteritems():
                savename='SNR_{0}'.format(k)
                store[savename]=v
        
        if self.sigma is not None:
            for k,v in self.sigma.iteritems():
                savename='sigma_{0}'.format(k)
                store[savename]=v
        
        store.close()
    
    def alt_resample(self, altrange, verbose=False):
        #takes a pandas dataframe generated by lnc and resamples on regular
        #intervals in altitude and resets the limits of the set
        #note: limits of altrange must be within original limits of altitude data
        
        if verbose:
            print 'Altitude step resampling in progress ...'
        

        self.BR=resample_cols(self.BR,altrange,verbose)
        self.PR=resample_cols(self.PR,altrange,verbose)
        #resample range corrected data
        if self.MSK is not None:
            self.MSK=resample_cols(self.MSK,altrange,verbose)
        else:
            if verbose:
                print "No Masked Profiles"
                
        if self.backscatter is not None:       
            self.backscatter=resample_cols(self.backscatter,altrange,verbose)
        else:
            if verbose:
                print "No Backscatter Profiles"

        if self.extinction is not None:       
            self.extinction=resample_cols(self.extinction,altrange,verbose)
        else:
            if verbose:
                print "No Extinction Profiles"

        if self.scenepanel is not None:
            paneldict={}
            for i in self.scenepanel.items:
                dftemp = self.scenepanel.loc[i]
                paneldict[i]=resample_cols(dftemp,altrange,verbose,method='ffill')                
            self.scenepanel=pan.Panel.from_dict(paneldict)
        else:
            if verbose:
                print "No Scene Analysis"
        
        if self.sigma is not None:
            self.calculate_sigma(winsize=10)
        
        if self.SNR is not None:
            self.calculate_SNR(winsize=10)
        print '... Done!'
        
        return self
    
    def time_resample(self, timestep=None, starttime=None,endtime=None, datamethod = 'mean',
                      sigma_winsize=10,SNR_winsize=10,verbose=False):
        #resamples a pandas dataframe generated by lnc_reader on a regular timestep
        #and optionally limits it to a preset time range
        #timestep must be in timeseries period format: numF where num=step size and
        #F = offset alias.  Ex: H = hours, M = minutes, S = seconds, L = millieconds
        if verbose:
                print 'Time step regularization in progress ...'
        
        dftemp=self.BR                  
        if starttime is not None:
            dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
        if endtime is not None:
            dftemp = dftemp.loc[dftemp.index<=endtime]  
        if timestep is not None:
            dftemp = dftemp.resample(timestep, how = datamethod)
        self.BR=dftemp

        dftemp=self.PR                  
        if starttime is not None:
            dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
        if endtime is not None:
            dftemp = dftemp.loc[dftemp.index<=endtime]  
        if timestep is not None:
            dftemp = dftemp.resample(timestep, how = datamethod)
        self.PR=dftemp
        
        if self.MSK is not None:
            dftemp=self.MSK                 
            if starttime is not None:
                dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
            if endtime is not None:
                dftemp = dftemp.loc[dftemp.index<=endtime]  
            if timestep is not None:
                dftemp = dftemp.resample(timestep, how = datamethod)
            self.MSK=dftemp
            
        if self.backscatter is not None:   
            dftemp=self.backscatter                 
            if starttime is not None:
                dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
            if endtime is not None:
                dftemp = dftemp.loc[dftemp.index<=endtime]  
            if timestep is not None:
                dftemp = dftemp.resample(timestep, how = datamethod)
            self.backscatter=dftemp
        
        if self.extinction is not None:
            dftemp=self.extinction                 
            if starttime is not None:
                dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
            if endtime is not None:
                dftemp = dftemp.loc[dftemp.index<=endtime]  
            if timestep is not None:
                dftemp = dftemp.resample(timestep, how = datamethod)
            self.extinction=dftemp
        
        if self.scenepanel is not None:
            paneltemp=self.scenepanel
            panelout=pan.Panel(items=paneltemp.items,major_axis=paneltemp.major_axis,
                               minor_axis=paneltemp.minor_axis)
            for i in paneltemp.items:
                dftemp=paneltemp[i]
                if starttime is not None:
                    dftemp = dftemp.loc[dftemp.index>=starttime]                                                     
                if endtime is not None:
                    dftemp = dftemp.loc[dftemp.index<=endtime]  
                if timestep is not None:
                    dftemp = dftemp.resample(timestep, how ='ffill')
                panelout[i]=dftemp
            self.scenepanel=panelout
        if verbose:
            print '... Done!'
        
        if self.sigma:
            self.calculate_sigma(winsize=sigma_winsize)
        
        if self.SNR:
            self.calculate_SNR(winsize=SNR_winsize)
                
        return self
        
    def calculate_sigma(self,winsize=10,verbose=False, datatypes=['all']):
        """
        Calculates stnadard deviations for LNC data
        
        inputs:
        num profs = number of vertical profiles to average together, defaults to 1
        datatypes = list of data types to callculate sigma for.  Could be 
                    'raw','rsq','NRB','depolrat', or 'all'
        
        output:
        self.sigma = a dict of pandas dataframes with datatype keys containing 
                    standard deviation values
        
        """
        
        if verbose:
            print "Calculating sigma"
        
        datasets=[]
        sigmadict={}
            
        for d in datatypes:
            if d=='BR' or d=='all':
                datasets.append(('BR',self.BR))
            if d=='PR' or d=='all':
                datasets.append(('PR',self.PR))
            if d=='backscatter' or d=='all':
                if self.backscatter is not None:
                    datasets.append(('backscatter',self.backscatter))
                elif verbose:
                    print "No Backscatter available for sigma calc"
            if d=='extinction' or d=='all':
                if self.extinction is not None:
                    datasets.append(('extinction',self.extinction))
                elif verbose:
                    print "No extinction available for sigma calc"
        
        for dset_name,dset in datasets: 
            sigmadict[dset_name] = []
            if verbose:
                print "Calculating sigma values for {0}".format(dset_name)
            tempdat=dset
            stdarray=pan.DataFrame(genfilt(tempdat,np.std,winsize),index=tempdat.index,
                                   columns=tempdat.columns)
            sigmadict[dset_name]=stdarray
                
        self.sigma=sigmadict

        if verbose:
            print "Sigma calculation done!"
            
        return self
    
    def calculate_SNR(self,winsize=10,verbose=False, datatypes=['all']):
        """
        Calculates signal to noise ratios for LNC data
        
        inputs:
        dfin = a pandas dataframe
        bg_alt = altitude above which signal is assumed to be purely background
                 if empty, topmost 100 data points are used
        num profs = number of vertical profiles to average together, defaults to 1
        datatypes = list of data types to callculate SNR for.  Could be 
                    'raw','rsq','NRB','depolrat', or 'all'
        
        output:
        self.SNR = a dict of pandas dataframes with datatype keys containing 
                    SNR values
        
        """
        
        if verbose:
            print "Calculating SNR"
        
        datasets=[]
        SNRdict={}
            
        for d in datatypes:
            if d=='BR' or d=='all':
                datasets.append(('BR',self.BR))
            if d=='PR' or d=='all':
                datasets.append(('PR',self.PR))
            if d=='backscatter' or d=='all':
                if self.backscatter is not None:
                    datasets.append(('backscatter',self.backscatter))
                elif verbose:
                    print "No Backscatter available for SNR calc"
            if d=='extinction' or d=='all':
                if self.extinction is not None:
                    datasets.append(('extinction',self.extinction))
                elif verbose:
                    print "No extinction available for SNR calc"
        
        for dset_name,dset in datasets: 
            stdarray=pan.DataFrame(genfilt(dset,np.std,winsize),index=dset.index,
                                   columns=dset.columns)
            meanarray=pan.DataFrame(genfilt(dset,np.mean,winsize),index=dset.index,
                                    columns=dset.columns)
            SNRtemp=(meanarray/stdarray).fillna(0.0)
            SNRdict[dset_name]=SNRtemp
            
        self.SNR=SNRdict

        if verbose:
            print "SNR calculation done!"
            
        return self
    
    def calc_all(self,winsize=10,verbose=False):
        """
        calculates all uncalculated fields for an LNC object
        """
        
        if self.sigma is None:
            self.calculate_sigma(winsize=winsize,verbose=verbose)
            
        if self.SNR is None:
            self.calculate_SNR(winsize=winsize,verbose=verbose)
        
        return self
   
        
def set_dir(titlestring):
     
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
     
    filenames = tkFileDialog.askopenfilename(title=titlestring, filetypes=[filetype],multiple='True')
    
    #do nothing if already a python list
    if filenames == "": 
        print "You didn't open anything!"  
        return
    
    root.destroy()
    
    if isinstance(filenames,list):
        result = filenames   
    elif isinstance(filenames,tuple): 
        result = list(filenames)
    else:
        #http://docs.python.org/library/re.html
        #the re should match: {text and white space in brackets} AND anynonwhitespacetokens
        #*? is a non-greedy match for any character sequence
        #\S is non white space
    
        #split filenames string up into a proper python list
        result = re.findall("{.*?}|\S+",filenames)
    
        #remove any {} characters from the start and end of the file names
        result = [ re.sub("^{|}$","",i) for i in result ] 
    result.sort()
    return result

def resample_cols(dfin,newcols,verbose=False,method='interp'):
    oldcols=dfin.columns
        
    mincol=oldcols[0]
    maxcol=oldcols[-1]
    
    if mincol>newcols[0]:
        newcols=sorted([c for c in newcols if c>=mincol])
        if verbose:
            print "WARNING: Minimum column value reset to {0}".format(newcols[0])
        
    if maxcol<newcols[-1]:
        newcols=sorted([c for c in newcols if c<=maxcol])
        if verbose:
            print "WARNING: Maximum column value reset to {0}".format(newcols[-1])
    
    if len(newcols)==len(oldcols) and all(newcols==oldcols):
        return dfin
    else:
        newvalues=[]   
        for row in dfin.iterrows():
            if method=='interp':
                f=interp1d(oldcols,row[1].values)
                newvalues.append(f(newcols))
            elif method=='ffill':
                newrow=[]
                for col in newcols:
                    edgeval=row[1].groupby(row[1].index<=col).groups[True][-1]
                    newrow.append(row[1][edgeval])
                newvalues.append(newrow)    
            elif method=='bfill':
                newrow=[]
                for col in newcols:
                    edgeval=row[1].groupby(row[1].index<=col).groups[False][0]
                    newrow.append(row[1][edgeval])
                newvalues.append(newrow)            
        dfout=pan.DataFrame(data=newvalues,index=dfin.index,columns=newcols)       
        return dfout
        
    


def partdepolratcalc(depolin,beta_parallel,beta_mol,moldepolrat=0.0035):
    
    #default moldepolrat: narrow double filter allows only Cabannes line (see SPIE proc reference)
#    A = moldepolrat/(1+moldepolrat)
#    partdepolrat=(depolin*beta_parallel-A*beta_mol)/(beta_parallel-A*beta_mol)
    beta_p=beta_parallel-beta_mol    
    partdepolrat=(beta_parallel*depolin - beta_mol*moldepolrat)/beta_p
    
    return partdepolrat

def buffered_array(data,(x,y)):

    #create buffer around dataframe
    datashape = np.shape(data)
    
    b = int(np.ceil(y/2))
    t = y-b
    
    if len(datashape)==1:
        rows=datashape[0]
        newsize=(rows+y)
        newarray=np.empty(newsize)
        (newrows)=newarray.shape
        newarray[b:-t]=data
        newarray[:b]=data[:b]
        newarray[-t:]=data[-t:]
        
    else:
        rows=datashape[0]
        columns=datashape[1]    
        #simply copy first and last values to fill in buffers
        if x > 1:
            l = int(np.ceil(x/2))
            r = x-l
            newsize=(rows+x,columns+y)
            newarray=np.empty(newsize)
            (newrows,newcolums)=newarray.shape
            newarray[:l,b:-t]=data[0,:]
            newarray[l:-r,b:-t]=data
            newarray[-r:,b:-t]=data[-1,:]
        else:
            l=0
            r=0
            newsize=(rows,columns+y)
            newarray=np.empty(newsize)
            (newrows,newcolums)=newarray.shape
            newarray[:,b:-t]=data
    
        newarray[:,:b]=newarray[:,b:2*b]
        newarray[:,-t:]=newarray[:,-2*t:-t] 
    
    return newarray

def SNR_mask_depol(LNCin,**kwargs):
    
    SNRthreshold=kwargs.get('SNRthreshold',3)
    numprofs=kwargs.get('numprofs',1)
    nopassval=kwargs.get('nopassval',float('nan'))
    inplace=kwargs.get('inplace',False)
    recalc=kwargs.get('recalc',False)
    
    if recalc or not LNCin.SNR:
        LNCin = LNCin.calculate_SNR(numprofs,datatypes='BR')
  
    #start by creating mask where areas that fall below SNRthreshold are zeroed out
    SNRmask = LNCin.SNR[datatype][0]>=SNRthreshold
    
    if inplace:
        LNCout=LNCin
    else:
        LNCout=deepcopy(LNCin)
           
    if LNCout.PR:
        LNCout.PR=LNCin.PR*SNRmask
        LNCout.PR.replace(0,nopassval,inplace=True)
    return LNCout 

def SNR_mask_scene(LNCin,**kwargs):
    
    SNRthreshold=kwargs.get('SNRthreshold',3)
    numprofs=kwargs.get('numprofs',1)
    bg_alt=kwargs.get('bg_alt',None)
    inplace=kwargs.get('inplace',False)
    recalc=kwargs.get('recalc',False)
    datatype=kwargs.get('datatype','NRB')
    nopassval=kwargs.get('nopassval',None)
    
    if nopassval is None:
        nopassval={'Base':np.nan,
                   'Top':np.nan,
                   'Lidar_Ratio':0.0,
                   'Delta':0.0,
                   'Type':'Insufficient Signal',
                   'Sub-Type':'Insufficient Signal',
                   'Depol':0.0,
                   'colormask':9}
    
    if recalc or LNCin.SNR is None:
        LNCin = LNCin.calculate_SNR(bg_alt,numprofs,datatypes=[datatype])
  
    #start by creating mask where areas that fall below SNRthreshold are zeroed out
    SNRmask = (LNCin.SNR[datatype][0]>=SNRthreshold)#|(LNCin.scenepanel[0]['Type']=='Clear Air')
#    SNRmask.replace(False,np.nan,inplace=True)
    
    if inplace:
        LNCout=LNCin
    else:
        LNCout=deepcopy(LNCin)
           
    if LNCout.scenepanel:
        for sceneitem,maskval in nopassval.iteritems():
            LNCout.scenepanel[0][sceneitem]=LNCin.scenepanel[0][sceneitem]*SNRmask
            LNCout.scenepanel[0][sceneitem].fillna(nopassval,inplace=True)
    return LNCout 
    
def SNR_mask_all(LNCin,**kwargs):
    
    SNRthreshold=kwargs.get('SNRthreshold',3)
    numprofs=kwargs.get('numprofs',1)
    bg_alt=kwargs.get('bg_alt',None)
    nopassval=kwargs.get('nopassval',float('nan'))
    inplace=kwargs.get('inplace',False)
    recalc=kwargs.get('recalc',False)
    datatype=kwargs.get('datatype','NRB')
    if recalc or not LNCin.SNR:
        LNCin = LNCin.calculate_SNR(bg_alt,numprofs,datatype=['data'])
  
    #start by creating mask where areas that fall below SNRthreshold are zeroed out
    SNRmask = LNCin.SNR[datatype][0]>=SNRthreshold
    
    if inplace:
        LNCout=LNCin
    else:
        LNCout=deepcopy(LNCin)
    
    LNCout.BR=LNCin.BR*SNRmask
    LNCout.BR.replace(0,nopassval,inplace=True)
    if LNCout.PR is not None:
        LNCout.PR=LNCin.PR*SNRmask
        LNCout.PR.replace(0,nopassval,inplace=True)        
    if LNCout.MSK is not None:
        LNCout.MSK=LNCin.MSK*SNRmask
        LNCout.MSK.replace(0,nopassval,inplace=True)      
    if LNCout.backscatter is not None:
        LNCout.backscatter=LNCin.backscatter*SNRmask
        LNCout.backscatter.replace(0,nopassval,inplace=True) 
    if LNCout.extinction is not None:
        LNCout.extinction=LNCin.extinction*SNRmask
        LNCout.extinction.replace(0,nopassval,inplace=True)   
    return LNCout  
    
def BR_mask_create(dfin,**kwargs):
    
    """
    generates threshold altitudes to avoids spurious results by removing all 
    data beyond strong signal spikes
    
    """
    BRthreshold=kwargs.get('BRthreshold',3)
    BRmin=kwargs.get('BRmin',0.05)
    minalt=kwargs.get('minalt',0.150)
    numprofs=kwargs.get('numprofs',1)
    winsize=kwargs.get('winsize',5)
    
    #start by creating array of threshold altitudes and masking BR copol
    #creates a new array with buffers to account for numprofs, winsize
    
    data = dfin.values
    altrange=dfin.columns.values
    (rows,columns) = data.shape
    minalt_index=np.where(altrange>=minalt)[0][0]
    newarray = buffered_array(data,(numprofs,winsize))
    (newrows,newcolums) = newarray.shape

    #set default values for cutoff to maximum altitude 
    threshalts=np.ones(len(dfin.index))*altrange[-1]
   
    for r in range(rows):
        tempprof=np.mean(newarray[r:r+numprofs],axis=0)
        for c in np.arange(minalt_index,columns):
            tempval = np.mean(tempprof[c:c+winsize])
            if tempval >= BRthreshold:
                for c2 in np.arange(c,columns):
                    tempval = np.mean(tempprof[c2:c2+winsize])
                    if tempval <= BRmin:
                        threshalts[r]=altrange[c2]
                        break 
    threshseries=pan.Series(data=threshalts,index=dfin.index)
    return threshseries

def BR_mask_apply(dfin,threshseries,nopassval=np.nan,inplace=True):
    
    if inplace:
        dfout=dfin
    else:
        dfout=deepcopy(dfin)
    altvals = dfin.columns.values    
    for r in dfin.index:
        tempval=[x for x in altvals if x>=threshseries.ix[r]][0]
        dfout.ix[r,tempval:]=nopassval
    return dfout

def BR_mask_all(LNCin,**kwargs):
    """
        uses a list of threshold altitudes, or generates one based on kwargs
        and applies it to all data sets within an LNC class object
    """
    
    threshseries=kwargs.get('threshseries',None)
    BRmasktype=kwargs.get('BRmasktype','profile')
    BRthreshold=kwargs.get('BRthreshold',3)
    BRmin=kwargs.get('BRmin',0.5)
    minalt=kwargs.get('minalt',0.150)
    numprofs=kwargs.get('numprofs',1)
    winsize=kwargs.get('winsize',3)
    nopassval=kwargs.get('nopassval',np.nan)
    inplace=kwargs.get('inplace',True)
    
    if inplace:
        LNCout=LNCin
    else:
        LNCout=deepcopy(LNCin)
    
    if threshseries is None:
        threshkwargs= {'BRthreshold':BRthreshold,'masktype':BRmasktype,'BRmin':BRmin,'minalt':minalt,
                       'numprofs':numprofs,'winsize':winsize,'nopassval':nopassval}
        threshseries=BR_mask_create(LNCout.BR,**threshkwargs)
    
    LNCout.BR=BR_mask_apply(LNCout.BR,threshseries)
    LNCout.PR=BR_mask_apply(LNCout.PR,threshseries)

    if LNCout.backscatter is not None:
        LNCout.backscatter=BR_mask_apply(LNCout.backscatter,threshseries)

    if LNCout.extinction is not None:
        LNCout.extinction=BR_mask_apply(LNCout.extinction,threshseries)
    
    if LNCout.scenepanel is not None:
        tempscene=LNCout.scenepanel
        tempscene['Type']=BR_mask_apply(tempscene['Type'],threshseries,nopassval='Insufficient Signal')
        tempscene['Sub-Type']=BR_mask_apply(tempscene['Sub-Type'],threshseries,nopassval='Insufficient Signal')
        tempscene['colormask']=BR_mask_apply(tempscene['colormask'],threshseries,nopassval=9)
        tempscene['Lidar_Ratio']=BR_mask_apply(tempscene['Lidar_Ratio'],threshseries,nopassval=np.nan)
        tempscene['Depol']=BR_mask_apply(tempscene['Depol'],threshseries,nopassval=np.nan)
        tempscene['Delta']=BR_mask_apply(tempscene['Delta'],threshseries,nopassval=np.nan)
        tempscene['Base']=BR_mask_apply(tempscene['Base'],threshseries,nopassval=np.nan)
        tempscene['Top']=BR_mask_apply(tempscene['Top'],threshseries,nopassval=np.nan)
    return LNCout
    


if __name__ == '__main__':

    olddir = os.getcwd()
    delta = 0.1

    os.chdir('K:\CORALNet\ASCII_Files\Smoke2012\UBC\July')
    
    BRfilename = get_files('Select BR pickle file',filetype=('.pickle','*.pickle'))[0]
    PRfilename = get_files('Select PR pickle file',filetype=('.pickle','*.pickle'))[0]
    MSKfilename = get_files('Select MSK pickle file',filetype=('.pickle','*.pickle'))[0]
    
    print 'Testing LNC class functions'
    
    print 'Import LNC data from .pickle file'
    
    LNCtest = LNC()
    LNCtest.frompickle(BRfilename,PRfilename,MSKfilename,altcor=True)
    
    print 'Done'
    
    print 'Calculate all corrections'
#    
    LNCtest.calc_all(verbose=False)
    
    LNCtest.save_to_HDF('UBC_20120703_PR532-UBC_20120720-proc_v3.h5')
#    os.chdir(olddir)

    os.chdir(olddir)

     
