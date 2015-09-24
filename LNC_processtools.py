# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 10:54:41 2014

@author: dashamstyr
"""
#from __future__ import absolute_import
import os,sys,site
home=os.environ['homepath']
site.addsitedir('{0}\\Dropbox\\Python_Scripts\\GIT_Repos\\'.format(home))

import pandas as pan
import numpy as np
from copy import deepcopy
from itertools import groupby
from scipy import signal
from matplotlib import pyplot as plt
import operator
import inversiontools as itools
from scipy import optimize as opt
from scipy.ndimage.filters import generic_filter as genfilt

import LNC_tools as ltools
import LNC_plot3 as lplot


def molecular_detect(LNCin,**kwargs):    
    
    wave=kwargs.get('wave',1064.0)
    winsize=kwargs.get('winsize',30) 
    deltathresh=kwargs.get('deltathresh',0.001)
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','testmolecular.h5')
    
        
    LNCin=LNCin.calc_all(winsize=winsize)  
    BRin=LNCin.BR    
    sigmain=LNCin.sigma['BR']
    
    #Step 1: calculate molecular profile
    z=BRin.columns.values
    altstep=z[1]-z[0]  #assumes regular altitude steps throughout
#    Pmol=itools.molprof(z,wave)
#    #LNC uses backscatter ratio, so this is the quantity used for comparison
#    Pmol_rat=Pmol['vals']/Pmol['beta_R']
    panelout=pan.Panel(major_axis=BRin.index,minor_axis=['Base','Top'])
    
    for proftime in BRin.index:    
    #Step 2:extract BR profile
        tempprof=BRin.ix[proftime]
        
        if np.isnan(tempprof).all():
            panelout.loc['Layer0',proftime]=[np.nan,np.nan]
            continue
        else:
        #Step 3: calculate ratio of profile to molecular and use smoothing window to reduce noise
#        temprat=Pmol_rat.div(tempprof)
#        coef=pan.Series(genfilt(temprat,np.mean,winsize),index=temprat.index)
#        #Step 4: Result from 3 is profile of multiplying factor.Use this to calculate variance of residuals
#        rawvariance=(1/(z**2.0)*(tempprof-(Pmol_rat/coef)))**2.0
#
#        variance=pan.Series(genfilt(rawvariance,np.mean,winsize),index=rawvariance.index)
#
#        sigmaprof=sigmain.ix[proftime]
        #Step 5: Regions where variance is below threshold multiple of noise varaince(sigma squared)
        #identified as molecular regions
#        tempmask=pan.Series(z,index=[v<=varthresh*s**2 for v,s in zip(variance,sigmaprof)])
            tempmask=pan.Series(z,index=tempprof.values)        
            tempgroups=tempmask.groupby(level=0)
            tempalts=[]
            for g in tempgroups:
                if abs(g[0]-1.0)<=deltathresh:
                    tempalts.extend(g[1].values)
            
            tempcounts=[int(round(x/altstep)) for x in tempalts]
            tempcounts.sort()                
            n=0
    #        oldlayer=[]
            for key,count in groupby(enumerate(tempcounts), lambda (i,x):i-x):                    
                layercount=map(operator.itemgetter(1),count)
    #                    if oldlayer and layercount[0]-oldlayer[-1]<=2:
    #                        layercount+=oldlayer
                layercount.sort()
    #            if len(layercount)==1:
    #                oldlayer+=layercount
    #            else:
                layeralt=[x*altstep for x in layercount]
    #                    if n!=0 and layeralt[0]<=panelout.loc['Layer{0}'.format(n-1),proftime][0]:
    #                        n=n-1
                panelname='Layer{0}'.format(n)
                panelout.loc[panelname,proftime]=[layeralt[0],layeralt[-1]]
                n+=1
#                    oldlayer=layercount
    if savefile:
        store=pan.HDFStore(savefilename)
        store['molecular']=panelout
        
    return panelout
        

def PBL_detect(LNCin,**kwargs):
    
    """
    """    
    wavelet=kwargs.get('wavelet',dog)
    widths=kwargs.get('widths',np.arange(3,10))
    layerwidth=kwargs.get('layerwidth',4)
    layer_min=kwargs.get('layer_min',None)
    mol_min=kwargs.get('mol_min',None)
    winsize=kwargs.get('winsize',10) 
    
    rawdata=LNCin.BR
    
    PBLout = pan.Series(index=rawdata.index)
    z=rawdata.columns.values
        
    for i in rawdata.index:
        tempprof=rawdata.ix[i] 
        if np.isnan(tempprof).all():
            PBLout.ix[i]=np.nan
            continue
        else:
            tempcwt=signal.cwt(tempprof,wavelet,widths)
            tempmin=maxmin(tempcwt,widths,np.less)
            
            #use only width of layerwidth for now            
            minloc=[minval[1] for minval in tempmin if minval[0]==layerwidth]
            
            #PBL height can only be negative (min) values and edge effects are removed
            #by remving edges layerwidth in size
            try:
                tempalt=z[minloc[1]]
            except IndexError:
                tempalt=z[minloc[0]]
            
            
#            layerwidthindex=np.where(widths==layerwidth)[0]
#            CWTvals=tempcwt[layerwidthindex,:][0]
#            CWTminvals=CWTvals[minloc]
#            CWTminalts=z[minloc]
#            
            tempmolht=mol_min.loc[i] 
            templayerht=layer_min.loc[i]
            
            edgealt=z[layerwidth-1]
            
            if np.isnan(templayerht) or templayerht<edgealt:
                templayeralt=z[-1]
            else:
                templayeralt=templayerht
            
            if np.isnan(tempmolht) or tempmolht<edgealt:
                tempmolalt=z[-1]
            else:
                tempmolalt=tempmolht
            
            maxalt=np.min([tempmolalt,templayeralt])
            
            if tempalt<=maxalt:
                PBLout.ix[i]=tempalt
            else:
                PBLout.ix[i]=maxalt
#            try:
#                PBLout.ix[i]=min([v for v in minalts if v<=maxalt])                
#            except ValueError:
#                PBLout.ix[i]=maxalt
                    
#            if PBLval is not None:
#                PBLout.ix[i]=CWTminalts[np.where(CWTminvals==PBLval)]
            
    return PBLout

def maxmin(arrayin,widths,f):
    #find all local maxima or minina from each row of an array and return array
    #of index values for each max or min
    
    arrayout=[]
    for n in range(len(arrayin)):
        temp=signal.argrelextrema(arrayin[n],f)[0]        
        for t in temp:
            arrayout.append((widths[n],t))
        
    return arrayout

def dog(points,a):
    y_out=[]
    for x in range(points):
        x=x-points/2
        y=(np.exp(-x**2/(2*a**2))*-x)/(np.sqrt(2*np.pi)*a**3)
        y_out.append(y)
    return y_out    


def find_layers(LNCin,**kwargs):
    """
    takes an LNC class object and process it, one profile at a time, to estimate
    bottom, peak,and top for each layer within the 2-D dataset
    
    inputs:
    LNCin = an MPL-class object to be proken into layers
    
    kwargs:
    wavelet=type of wavelet to use to find layer edges.  default:signal.ricker
    widths=Range of wavelet widths to feed into CWT.  default:[2]
    layerwidth=wavelet width to use for final layer ID. default:2
    bg_alt=altitude to mark as particulate-free for background calc. defualt:[]
    noisethresh=threshold level to mark a layer above background noise. default:3
    cloudthresh=threshold signal level to mark a layer as cloud for (water,ice). default: (1.0,0.4)
    datatype=type of profile in MPL object to use for layers.  default:'data'
    savefile=boolean for whetehr to save results. default:False
    savefilename=name to save file under.  default:'testlayers.h5'
    
    
    Outputs:
    panelout = a pandas panel object with three axes:
        major-axis: datetime of individual profiles
        minor-axis: layer info ['Base','Peak','Top','Delta','Depol','Type',
                                'Sub-Type','Lidar_Ratio']
        columns: Layer number (e.g. 'Layer1')
    """
    
    #if MPLin does not have all necessary processed data,generate it
    wavelet=kwargs.get('wavelet',signal.ricker)
    widths=kwargs.get('widths',np.arange(2,5))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    noisethresh=kwargs.get('noisethresh',3)
    minwidth=kwargs.get('minwidth',4)
    minsep=kwargs.get('minsep',2)
    sigma0=kwargs.get('sigma0',None)
    depolsigma0=kwargs.get('depolsigma0',None)
    cloudthresh=kwargs.get('cloudthresh',(10.0,2.0))
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    datatype=kwargs.get('datatype','NRB')
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','testlayers.h5')    
    winsize=kwargs.get('winsize',10) 
    bg_alt=kwargs.get('bg_alt',None)
    udefrat=kwargs.get('udefrat',0)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    
    if bg_alt is None:
        bg_alt=LNCin.BR.columns.values[-10]
    
    LNCin=LNCin.calc_all(winsize=winsize) 
    #use raw data from co-polarized channel (not r-squared corrected) to find layers
    rawdata=LNCin.BR
    rawdepol=LNCin.PR
    rawsigma=LNCin.sigma['BR']
    depolsigma=LNCin.sigma['PR']
    panelout=pan.Panel(major_axis=rawdata.index,minor_axis=['Base','Peak','Top',
    'Delta','Depol','Type','Sub-Type','Lidar_Ratio'])
    
    for i in rawdata.index:
        tempprof=rawdata.ix[i]
        
        if np.isnan(tempprof).all():
            panelname='Layer0'
            panelout.loc[panelname,i,'Base']=np.nan
            panelout.loc[panelname,i,'Peak']=np.nan
            panelout.loc[panelname,i,'Top']=np.nan
            panelout.loc[panelname,i,'Delta']=np.nan
            panelout.loc[panelname,i,'Depol']=np.nan
            panelout.loc[panelname,i,'Type']=np.nan
            panelout.loc[panelname,i,'Sub-Type']=np.nan
            panelout.loc[panelname,i,'Lidar_Ratio']=np.nan
            continue
        else:            
            tempdepolprof=rawdepol.ix[i]
            tempsigma=rawsigma.ix[i]
            tempdepolsigma=depolsigma.ix[i]
    #        tempdepolratprof=rawdepolrat.ix[i]
    #        z=tempprof.index
            #set baseline noise level based on 
            if sigma0:
                tempsigma0=sigma0
            else:
                tempsigma0=np.mean(tempsigma.ix[bg_alt:])
            
            if depolsigma0:
                tempdepolsigma0=depolsigma0
            else:
                tempdepolsigma0=np.mean(tempdepolsigma.ix[bg_alt:])
            
            temp_cwt=signal.cwt(tempprof,wavelet,widths)
            tempmax=maxmin(temp_cwt,widths,np.greater)
            tempmin=maxmin(temp_cwt,widths,np.less)
            
            #use only width of 2 for now
            
            minloc=[minval[1] for minval in tempmin if minval[0]==CWTwidth]
            maxloc=[maxval[1] for maxval in tempmax if maxval[0]==CWTwidth]
            filterkwargs={'thresh':noisethresh,'minwidth':minwidth,'minsep':minsep,
                          'depolwidths':widths,'depollayerwidth':minwidth,
                          'depolwavelet':wavelet}
            templayers=layer_filter(tempprof,tempsigma0,tempdepolprof,tempdepolsigma0,
                                    maxloc,minloc,datatype,**filterkwargs)
            
            for n in range(len(templayers)):
                indices=templayers[n]
                
                minalt=indices[0]
                peakalt=indices[1]
                maxalt=indices[2]
                delta=indices[3]
                meandepolrat=indices[4]
                panelname='Layer{0}'.format(n)
    #            layerdepolratprof=tempdepolratprof.ix[minalt:maxalt]
    #            meandepolrat=np.mean(layerdepolratprof)
                peakval=tempprof.ix[peakalt]
                if peakval <= 0.0 or meandepolrat <= 0.0:
                    layertype='Insufficient Signal'
                    layersubtype='Insufficient Signal'
                    layerratio=udefrat
                elif meandepolrat >= 0.5:
                    layertype='Cloud'  
                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                           icethresh=icethresh)    
                elif minalt>= maxaeroalt:
                    layertype='Cloud'  
                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                           icethresh=icethresh)
                elif peakval >= cloudthresh[0]:
                    layertype='Cloud'  
                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                           icethresh=icethresh)    
                elif peakval >= cloudthresh[1] and meandepolrat>np.mean([waterthresh,icethresh]):
                    layertype='Cloud'  
                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                           icethresh=icethresh)
                else:
                    layertype='Aerosol'
                    layersubtype,layerratio=aerosoltypefilter(meandepolrat,smokethresh=smokethresh,
                                                              dustthresh=dustthresh)
                                    
                panelout.loc[panelname,i,'Base']=minalt
                panelout.loc[panelname,i,'Peak']=peakalt
                panelout.loc[panelname,i,'Top']=maxalt
                panelout.loc[panelname,i,'Delta']=delta
                panelout.loc[panelname,i,'Depol']=meandepolrat
                panelout.loc[panelname,i,'Type']=layertype
                panelout.loc[panelname,i,'Sub-Type']=layersubtype
                panelout.loc[panelname,i,'Lidar_Ratio']=layerratio
    
    if savefile:
        store=pan.HDFStore(savefilename)
        store['layers']=panelout
        
    return panelout

def layer_filter(prof,sigma0,depolratprof,depolratsigma0,maxiloc,miniloc,datatype,**kwargs):
    """
    takes a profile and a list of local maxima and minima from CWT analysis and calculates
    layer edges and peaks while filtering out peaks for which the delta from 
    edge to peak is less than some multiple of the shot noise from background and
    dark current
    
    once a layer is defined, it is then investigated for variations in depol ratio
    If significant variations exist, the layer is firther divided into sub-layers 
    based on these results
    
    inputs:
    prof - a pandas series represeting a single profile of lidar returns with altitude
    depolprof - a pandas series representing a profile of depol ratios with altitude
    maxiloc - a list of maximum index values from the CWT results at a given wavelet width
            represent the peaks of a given layer
    miniloc - a list of minimum index values from the CWT results at a given wavelet width
            represent the edges of a given layer
    sigma0 - baseline noise level for the profile, if empty it is calculated
    thresh - difference between peak and edge of a layer must exceeed this 
             multiple of sigma0 to be counted.  default: 3
    
    """
    thresh=kwargs.get('thresh',3)
    minwidth=kwargs.get('minwidth',4)
    minsep=kwargs.get('minsep',2)
    depolwidths=kwargs.get('depolwidths',np.arange(2,5))
    depollayerwidth=kwargs.get('depollayerwidth',4)
    depolwavelet=kwargs.get('depolwavelet',signal.ricker)

    #Step 1: Calculate profile values at each edge and peak
    layers=[]
    n=0
    nextminloc=0
    
    altstep=prof.index[1]-prof.index[0]
    
    while n < len(maxiloc)-1:
        n+=1
        #note: due to buffering to obtain averages, first and last peaks are discarded
        peakloc=maxiloc[n]
        peakval=prof.iloc[peakloc]
        threshval=thresh*sigma0
        
        if peakval >= threshval:
            edge_below_list=[v for v in miniloc[nextminloc:] if v<peakloc]
            edge_above_list=[v for v in miniloc[nextminloc:] if v>peakloc]
            
            if not edge_above_list or not edge_below_list:
                continue
            #Step 3: Calculate delta signal between peak and lower edge (directly before)
            for edge_below in edge_below_list[::-1]:            
                delta_lower=prof.iloc[peakloc]-prof.iloc[edge_below]        
                #Step 4: Filter out false layers for which delta < thresh*signam0
                if delta_lower>threshval and edge_below>=1.0:
                    templowedge=edge_below
                    break
                else:
                    templowedge=None
                    #try to find upper edge where delta_upper exceeds threshold
            if templowedge is not None:
                for edge_above in edge_above_list:
                    delta_upper=prof.iloc[peakloc]-prof.iloc[edge_above]
                    if delta_upper>threshval and edge_above>=1.0:
                        #if upper edge is found, add indices of (lower,center,upper, maximum delta) to layers
    #                    temppeakval=np.max(prof.iloc[templowedge:edge_above])
    #                    temppeakloc=np.where(prof.values==temppeakval)
                        tempprof=prof.iloc[templowedge:edge_above]
                        if len(tempprof)>=minwidth:
                            tempdepolratprof=depolratprof.iloc[templowedge:edge_above]
        #                    delta=max(delta_lower,delta_upper)
                            depolkwargs={'widths':depolwidths,'layerwidth':depollayerwidth,
                                         'wavelet':depolwavelet}
                            depol_layers=find_depollayers(tempprof,tempdepolratprof,
                                                          depolratsigma0,**depolkwargs)
                            
                            layers+=depol_layers
    #                    layers.append((templowedge,temppeakloc,edge_above,max(delta_lower,delta_upper)))
                        try:
                            nextpeak=[p for p in maxiloc if p >edge_above][0]
                        except IndexError:
                            break
                        nextminloc=miniloc.index(edge_above)
                        n=maxiloc.index(nextpeak)-1
                        break
    return layers

def find_depollayers(copolprof,depolratprof,depolratsigma0,**kwargs):
    widths=kwargs.get('widths',np.arange(2,5))
    layerwidth=kwargs.get('layerwidth',2)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',1)
    
#    prevlayernum=kwargs.get('prevlayernum',0)
#    depolratSNR=calc_SNR_depolrat(depolratprof,depolratsigma=depolratsigma,signal0=signalsigma0)
    
    temp_cwt=signal.cwt(depolratprof,wavelet,widths)
    tempmax=maxmin(temp_cwt,widths,np.greater)
    tempmin=maxmin(temp_cwt,widths,np.less)
    z=depolratprof.index
    
    minloc=[minval[1] for minval in tempmin if minval[0]==layerwidth]
    maxloc=[maxval[1] for maxval in tempmax if maxval[0]==layerwidth]
    
    edgeloc=np.sort(minloc+maxloc)
    edgealts=[z[v] for v in edgeloc]
    templayers=depol_filter(depolratprof,copolprof,edgealts,depolratsigma0,noisethresh)
    
    return templayers
    
def depol_filter(depolprof,signalprof,edgealts,sigma0,thresh=1):
    """
    
    """
    #Step 2: Calculate profile values at each edge and peak
    
    base=depolprof.index[0]
    n=0
    layers=[]
    while n<len(edgealts):
        edge=edgealts[n]
        try:
            nextedge=edgealts[n+1]
        except IndexError:
            nextedge=depolprof.index[-1]
           
        meanbelow=np.mean(depolprof.ix[base:edge])
        meanabove=np.mean(depolprof.ix[edge:nextedge])
        depoldelta=abs(meanabove-meanbelow)
        if depoldelta>=thresh*sigma0:
            if n==len(edgealts):
                bottom=edge
                top=nextedge
                meandepol=meanabove
            else:
                bottom=base
                top=edge
                meandepol=meanbelow
            temppeakval=np.max(signalprof.ix[bottom:top])
            temppeakloc=np.where(signalprof.values==temppeakval)
            temppeakalt=signalprof.index[temppeakloc].values[0]
            delta_below=signalprof.ix[temppeakalt]-signalprof.ix[bottom]
            delta_above=signalprof.ix[temppeakalt]-signalprof.ix[top]
            signaldelta=max(delta_below,delta_above)
            layers.append([bottom,temppeakalt,top,signaldelta,meandepol])
            base=edge
            n+=1
        else:
            n+=1
    if not layers:
        bottom=depolprof.index[0]
        top=depolprof.index[-1]
        meandepol=np.mean(depolprof)
        temppeakval=np.max(signalprof)
        temppeakloc=np.where(signalprof.values==temppeakval)
        temppeakalt=signalprof.index[temppeakloc].values[0]
        delta_below=signalprof.ix[temppeakalt]-signalprof.ix[bottom]
        delta_above=signalprof.ix[temppeakalt]-signalprof.ix[top]
        signaldelta=max(delta_below,delta_above)
        layers.append([bottom,temppeakalt,top,signaldelta,meandepol])
           
    return layers

def layerprofplot(profin,layersin,numlayer=30):
    
    z=profin.index
    vals=profin.values
    numfigs=len(plt.get_fignums())
    fig=plt.figure(numfigs+1)
    ax=fig.add_subplot(111)
    ax.plot(vals,z)
    
    mcolors=['blue','red','green','yellow','orange','purple']
    if numlayer>len(layersin.ix['Base']):
        numlayer=len(layersin.ix['Base'])
        
    for n in range(numlayer):
        if n>(len(mcolors)-1):
            color=mcolors[int(np.floor(n/len(mcolors)))-1]
        else:
            color=mcolors[n]   
        templayer=layersin.iloc[:,n]
        if not np.isnan(templayer['Base']):
            ax.scatter(profin.ix[templayer.iloc[0]],templayer.iloc[0],c=color,marker='o')
            ax.scatter(profin.ix[templayer.iloc[1]],templayer.iloc[1],c=color,marker='x')
            ax.scatter(profin.ix[templayer.iloc[2]],templayer.iloc[2],c=color,marker='v')
            
    fig.canvas.draw()


def icewaterfilter(depolrat,**kwargs):
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    
    if depolrat <= waterthresh:
        typeout='Water Cloud'
        ratout=15.3
    elif depolrat <= icethresh:
        typeout='Mixed Cloud'
        ratout= 25.0
    else:
        typeout='Ice Cloud'
        ratout=50.0
    
    return typeout,ratout


def aerosoltypefilter(depolrat,**kwargs):
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    
    if depolrat <= smokethresh:
        typeout='Smoke / Urban'
        ratout=65.0
    elif depolrat <= dustthresh:
        typeout='Polluted Dust'
        ratout=60.0
    else:
        typeout='Dust'
        ratout=55.0
    
    return typeout,ratout

    
def colormask_fromdict(LNCin,pblin,molin,layersin):
    
    alts=LNCin.BR.columns
    times=LNCin.BR.index
    mask=pan.DataFrame(index=times,columns=alts)
    
    colordict={'Clear Air':0,
               'PBL':1,
               'Ice Cloud':2,
               'Water Cloud':3,
               'Mixed Cloud':4,
               'Dust':5,
               'Polluted Dust':6,
               'Smoke / Urban':7,
               'Unidentified Aerosol':8,
               'Insufficient Signal':9}
    
    
    for t in times:
        tempprof=pan.Series(index=alts)
        
        for m in molin.items:
            tempmol=molin.ix[m,t]
            tempminalt=tempmol.ix['Base']
            tempmaxalt=tempmol.ix['Top']
            tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=colordict['Clear Air']

        
        for l in layersin.items:
            templayer=layersin.ix[l,t]
            tempminalt=templayer.ix['Base']
            tempmaxalt=templayer.ix['Top']
            temptype=templayer.ix['Type']
            tempsubtype=templayer.ix['Sub-Type']
            try:
                tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=colordict[tempsubtype]
            except KeyError:
                tempprof[(tempprof.index>=tempminalt) & (tempprof.index<=tempmaxalt)]=8
                
        PBLht=pblin.ix[t]
        tempprof[tempprof.index<=PBLht]=1

        tempprof.fillna(value=8,inplace=True)
        mask.ix[t]=tempprof
    
    return mask,colordict


    
def findalllayers(**kwargs):
    LNCin=kwargs.get('LNCin',None)
    filename=kwargs.get('filename',None)
    timestep=kwargs.get('timestep','240S')
    bg_alt=kwargs.get('bg_alt',None)
    molthresh=kwargs.get('molthresh',1)
    winsize=kwargs.get('winsize',5)
    wavelet=kwargs.get('wavelet',signal.ricker)
    noisethresh=kwargs.get('noisethresh',0.4)
    cloudthresh=kwargs.get('cloudthresh',(10.0,2.0))
    CWTwidth=kwargs.get('CWTwidth',2)
    minwidth=kwargs.get('minwidth',4)
    layerCWTrange=kwargs.get('layerCWTrange',np.arange(2,5))
    PBLwavelet=kwargs.get('PBLwavelet',dog)
    PBLCWTrange=kwargs.get('PBLCWTrange',np.arange(2,10))
    PBLwidth=kwargs.get('PBLwidth',7)
    savemasks=kwargs.get('savemasks',False)
    savemaskname=kwargs.get('savemaskname','testmasksall.h5')
    sigma0=kwargs.get('sigma0',None)
    depolsigma0=kwargs.get('depolsigma0',None)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    
    if LNCin is None:
        if filename is not None:
            LNCin = ltools.LNC()
            LNCin.fromHDF(filename)
        else:
            print "No LNC file or filename provided"
            return
    
    LNCin.time_resample(timestep=timestep)
    LNCin.calc_all()     
    
    molecular=molecular_detect(LNCin,varthresh=molthresh, winsize=winsize)
    
    layerkwargs = {'wavelet':wavelet,'noisethresh':noisethresh,
                   'cloudthresh':cloudthresh,'CWTwidth':CWTwidth,
                   'widths':layerCWTrange,'minwidth':minwidth,'bg_alt':bg_alt,
                   'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                   'dustthresh':dustthresh,'sigma0':sigma0,'depolsigma0':depolsigma0,'maxaeroalt':maxaeroalt}
    layers=find_layers(LNCin,**layerkwargs)  
    mol_min=molecular.loc['Layer0']['Base']
    try:
        layer0_base=layers.loc['Layer0']['Base']
        layer0_top=layers.loc['Layer0']['Top']
        profmin=LNCin.BR.columns[-1]
        layer_min = pan.Series(data=[b if b>profmin else t for b,t in zip(layer0_base,layer0_top)],index=LNCin.BR.index)

    except KeyError:
        layer_min=LNCin.BR.columns[-1]
    
    if doPBL:    
        PBLkwargs = {'wavelet':PBLwavelet,'mol_min':mol_min,'layer_min':layer_min,
                     'widths':PBLCWTrange,'layerwidth':PBLwidth,'bg_alt':bg_alt}
        pbl=PBL_detect(LNCin,**PBLkwargs)
    else:
        pbl=pan.Series(0.0,index=LNCin.BR.index)
    
    if savemasks:    
        store=pan.HDFStore(savemaskname)
        store['molecular']=molecular
        store['layers']=layers
        store['PBL']=pbl
        store.close()
    
    dictout={'LNC':LNCin,'molecular':molecular,'layers':layers,'pbl':pbl}
    return dictout   

def scenemaker(layerdict,**kwargs):
    """
    Takes outputs from layer masking subroutines and generates a single pandas panel with all information
    needed to enter extinction processing phase
    
    Inputs:
    layerdict - dict objct containing the following key-value pairs
    
    kwargs:
    PBLrat - 
    
    """
    LNCin=layerdict['LNC']
    molecular=layerdict['molecular']
    layers=layerdict['layers']
    PBL=layerdict['pbl']
    
    PBLrat=kwargs.get('PBLrat',30.0)
    molrat=kwargs.get('molrat',0.0)
    moldepol=kwargs.get('moldepol',0.0035)
    udefrat=kwargs.get('udefrat',molrat)
    udefdepol=kwargs.get('udefdepol',moldepol)
    savefile=kwargs.get('savefile',False)
    savefilename=kwargs.get('savefilename','test.h5')
    colordict=kwargs.get('colordict',None)
    waterthresh=kwargs.get('waterthresh',0.10)
    icethresh=kwargs.get('icethresh',0.35)
    smokethresh=kwargs.get('smokethresh',0.10)
    dustthresh=kwargs.get('dustthresh',0.20)
    maxaeroalt=kwargs.get('maxaeroalt',10.0)
    cloudthresh=kwargs.get('cloudthresh',(10.0,2.0))
    minlayerwidth=kwargs.get('minlayerwidth',4)
    
    if colordict is None:
        colordict={'Clear Air':0,
                   'PBL':1,
                   'Ice Cloud':2,
                   'Water Cloud':3,
                   'Mixed Cloud':4,
                   'Dust':5,
                   'Polluted Dust':6,
                   'Smoke / Urban':7,
                   'Unidentified Aerosol':8,
                   'Insufficient Signal':9}
                   
    alts=LNCin.BR.columns
    times=LNCin.BR.index
    #initialize dataframes 
    typemask=pan.DataFrame(index=times,columns=alts)
    subtypemask=pan.DataFrame(index=times,columns=alts)
    lrat=pan.DataFrame(index=times,columns=alts)
    depolrat=pan.DataFrame(index=times,columns=alts)
    delta=pan.DataFrame(index=times,columns=alts)
    layerbase=pan.DataFrame(index=times,columns=alts)
    layertop=pan.DataFrame(index=times,columns=alts)
    colormask=pan.DataFrame(index=times,columns=alts)
    

    
    #fill dataframes with values defined by layerdict
    for t in times:
        if np.isnan(LNCin.BR.ix[t]).all():
            typemask.loc[t,:]=np.nan
            subtypemask.loc[t,:]=np.nan           
            lrat.loc[t,:]=np.nan
            depolrat.loc[t,:]=np.nan
            delta.loc[t,:]=np.nan
            layerbase.loc[t,:]=np.nan
            layertop.loc[t,:]=np.nan  
            colormask.loc[t,:]=np.nan
            continue
        else:
            tempmol=molecular.ix[:,t].dropna(axis=1,how='all')
            templayers=layers.ix[:,t].dropna(axis=1,how='all')
    
            #then assign molecular props to altitudes identified as such
            for m in tempmol.iteritems():
                base=m[1]['Base']
                top=m[1]['Top']
                typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]='Clear Air'
                subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]='Clear Air'            
                lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=molrat
                depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=moldepol
                delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=0.0
                layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
                layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top  
                colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict['Clear Air']
            #finally assign layer properties to each layer one by one
            for l in templayers.iteritems():
                base=l[1]['Base']
                top=l[1]['Top']
                temprat=l[1]['Lidar_Ratio']
                tempdelta=l[1]['Delta']
                temptype=l[1]['Type']
                tempsubtype=l[1]['Sub-Type']
                tempdepol=l[1]['Depol']
                typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]=temptype
                subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]=tempsubtype
                lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=temprat
                depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=tempdepol
                delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=tempdelta
                layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
                layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top
                colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict[tempsubtype]
            
            tempPBL=PBL.ix[t]
            #finally assign PBL props to altitudes below PBL height
            typemask.loc[t,typemask.columns<=tempPBL]='PBL'
            subtypemask.loc[t,subtypemask.columns<=tempPBL]='PBL'
            lrat.loc[t,lrat.columns<=tempPBL]=PBLrat
            depoltemp=LNCin.PR.ix[t]
            depolrat.loc[t,depolrat.columns<=tempPBL]=np.mean(depoltemp[depoltemp.index<tempPBL])
            layerbase.loc[t,layerbase.columns<=tempPBL]=alts[0]
            layertop.loc[t,layertop.columns<=tempPBL]=tempPBL
            colormask.loc[t,colormask.columns<=tempPBL]=colordict['PBL']
            
            #find and classify unidentified areas
            tempprof=typemask.loc[t]
            tempprof.fillna('Unidentified',inplace=True)
            tempmask=pan.Series(alts,index=[v=='Unidentified' for v in tempprof])
            tempgroups=tempmask.groupby(level=0) 
            #assuming uniform altitude steps
            altstep=alts[1]-alts[0]
            for g in tempgroups:
                if g[0]:
                    tempalts=g[1]
                    tempcounts= [int(round((x-tempalts.iloc[0])/altstep)) for x in tempalts]
                    for key,count in groupby(enumerate(tempcounts), lambda (i,x):i-x):                    
                        layercount=map(operator.itemgetter(1),count)
                        layeralt=[x*altstep+tempalts.iloc[0] for x in layercount]
                        if len(layeralt)==1:
                            continue
                        else:
                            base=layeralt[0]
                            top=layeralt[-1]
                            unidentifiedBR=LNCin.BR.loc[t,base:top]
                            unidentifiedPR=LNCin.PR.loc[t,base:top]
                            if len(unidentifiedBR) <= minlayerwidth:
                                continue
                            else:
                                peakval=max(unidentifiedBR)
                                meandepolrat=np.mean(unidentifiedPR)
                                layerdelta=np.max([peakval-unidentifiedBR.iloc[0],peakval-unidentifiedBR.iloc[-1]])
                                
                                if peakval <= 0.0 or meandepolrat <= 0.0:
                                    layertype='Insufficient Signal'
                                    layersubtype='Insufficient Signal'
                                    layerratio=udefrat
                                    meandepolrat=udefdepol
                                    
                                elif meandepolrat >= 0.5:
                                    layertype='Cloud'  
                                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                                           icethresh=icethresh)    
                                elif base>= maxaeroalt:
                                    layertype='Cloud'  
                                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                                           icethresh=icethresh)
                                elif peakval >= cloudthresh[0]:
                                    layertype='Cloud'  
                                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                                           icethresh=icethresh)    
                                elif peakval >= cloudthresh[1] and meandepolrat>np.mean([waterthresh,icethresh]):
                                    layertype='Cloud'  
                                    layersubtype,layerratio=icewaterfilter(meandepolrat,waterthresh=waterthresh,
                                                                           icethresh=icethresh)
                                else:
                                    layertype='Aerosol'
                                    layersubtype,layerratio=aerosoltypefilter(meandepolrat,smokethresh=smokethresh,
                                                                              dustthresh=dustthresh)
                                                                              
                                typemask.loc[t,(typemask.columns>=base)&(typemask.columns<=top)]=layertype
                                subtypemask.loc[t,(subtypemask.columns>=base)&(subtypemask.columns<=top)]=layersubtype
                                lrat.loc[t,(lrat.columns>=base)&(lrat.columns<=top)]=layerratio
                                depolrat.loc[t,(depolrat.columns>=base)&(depolrat.columns<=top)]=meandepolrat
                                delta.loc[t,(delta.columns>=base)&(delta.columns<=top)]=layerdelta
                                layerbase.loc[t,(layerbase.columns>=base)&(layerbase.columns<=top)]=base
                                layertop.loc[t,(layertop.columns>=base)&(layertop.columns<=top)]=top
                                colormask.loc[t,(colormask.columns>=base)&(colormask.columns<=top)]=colordict[layersubtype]
                            
    typemask.fillna(method='ffill',axis=1,inplace=True)
    subtypemask.fillna(method='ffill',axis=1,inplace=True) 
    lrat.fillna(method='ffill',axis=1,inplace=True)
    depolrat.fillna(method='ffill',axis=1,inplace=True)
    layerbase.fillna(method='ffill',axis=1,inplace=True)
    layertop.fillna(method='ffill',axis=1,inplace=True)
    colormask.fillna(method='ffill',axis=1,inplace=True)
    paneldict={'Type':typemask,'Sub-Type':subtypemask,'Lidar_Ratio':lrat,'Depol':depolrat,
               'Delta':delta,'Base':layerbase,'Top':layertop,'colormask':colormask}
    panelout=pan.Panel.from_dict(paneldict)
    
    LNCin.scenepanel=panelout

    if savefile:
        LNCin.save_to_HDF(savefilename)        

    return LNCin

def Tfromextinction(extinctprof): 
    alts=extinctprof.index
    tau=0.0
    oldz=alts[0]
    oldex=extinctprof.ix[oldz]
    for z in alts[1:]:
        deltaz=z-oldz
        tempex=(oldex+extinctprof.ix[z])/2.0
        tau+=-2*tempex*deltaz
        oldz=z
        oldex=extinctprof.ix[z]    
    Tout=np.exp(tau)
    return Tout
    
def Tfromprofile(prof,refprof=None,**kwargs):
    #requires profile containing at least [numvals] points on top and bottom
    #of molecular scattering free of aerosols
    wave=kwargs.get('wave',532.0)
    average=kwargs.get('average',True)
    numvals=kwargs.get('numvals',10)
    
    alts=prof.index.values
    if refprof is None:
        refprof=itools.molprof(alts,wave=wave)['vals']*alts**2
    
    if average:
        profabove=np.mean(prof.iloc[-numvals:])
        profbelow=np.mean(prof.iloc[:numvals])
        refabove=np.mean(refprof.iloc[-numvals:])
        refbelow=np.mean(refprof.iloc[:numvals])
    else:
        profabove=prof.iloc[-1]
        profbelow=prof.iloc[0]
        refabove=prof.iloc[-1]
        refbelow=prof.iloc[0]
    
    normrat=refbelow/profbelow
    Tout=profabove*normrat/refabove
    return Tout

def Tmolecular(minalt,maxalt,refprof=None,wave=532.0,altstep=0.030):
    alts=np.arange(minalt,maxalt+altstep,altstep)
    if refprof is None:
        refprof=itools.molprof(alts,wave=wave)['sigma_R']
    tau=0.0
    oldz=alts[0]
    oldex=refprof.ix[oldz]
    for z in alts[1:]:
        deltaz=z-oldz
        tempex=(oldex+refprof.ix[z])/2.0
        tau+=-2*tempex*deltaz
        oldz=z
        oldex=refprof.ix[z]    
    Tout=np.exp(tau)
    return Tout

def Tfrombounds(molbelow,molabove,wave=532.0, altstep=0.030):
    refalts=np.arange(molbelow.index[0],molabove.index[-1]+altstep,altstep)
    refprof=itools.molprof(refalts,wave=wave)['vals']*refalts**2
    refbelow=refprof.ix[molbelow.index]
    refabove=refprof.ix[molabove.index]        
    normrat=np.mean(refbelow.values)/np.mean(molbelow.values)
    Tout=np.mean(molabove.values)*normrat/np.mean(refabove.values)
    return Tout

def coefprofs(profin,lratprof,**kwargs):
    method=kwargs.get('method','klett2')
    refalt=kwargs.get('refalt',profin.index[-1])
    energy=kwargs.get('energy',1.0)
    calrange=kwargs.get('calrange',None)
    wave=kwargs.get('wave',2064.0)
    
    if pan.isnull(profin).all():
        backprof=pan.Series(0,index=profin.index)
        extprof=backprof
    else:
        gtemp=pan.groupby(profin,pan.isnull(profin))
        limalt=gtemp.groups[False][-1]
        if limalt<refalt:
            refalt=limalt
        if method=='fernald':
            if calrange is None:
                mincalalt=profin.index.values[:refalt][-20]
                calrange=[mincalalt,refalt]
            lratin=np.mean(lratprof)
            backprof,extprof=itools.fernald(profin,lratin,E=energy,calrange=calrange)
        elif method=='klett':
            backprof,extprof=itools.klett(profin,lratprof,r_m=refalt)
        elif method=='klett2':
            backprof,extprof=itools.klett2(profin,lratprof,r_m=refalt)
    return backprof,extprof
    
def basiccorrection(LNCin,**kwargs):    
    wave=kwargs.get('wave',1064.0)
    refalt=kwargs.get('refalt',None)
    calrange=kwargs.get('calrange',None)
    method=kwargs.get('method','klett2')
    lrat=kwargs.get('lrat',None)
        
    BRdat=LNCin.BR
    scenepanel=LNCin.scenepanel
    
    times=BRdat.index
    alts=BRdat.columns
    tempbackscatter=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    tempextinction=pan.DataFrame(data=np.nan,index=times,columns=alts,dtype='float')
    
    for t in times:
        tempprof=BRdat.ix[t]
        if np.isnan(tempprof).all():
            tempbackscatter.ix[t]=np.nan
            tempextinction.ix[t]=np.nan
            continue
        else:
            if refalt is None:
                refalt=tempprof.index.values[-1]
                gtemp=pan.groupby(tempprof,pan.isnull(tempprof))
                limalt=gtemp.groups[False][-1]
                if limalt<refalt:
                    refalt=limalt
            if lrat:
                templrat=pan.Series(data=lrat,index=alts)
            else:
                templrat=scenepanel.loc['Lidar_Ratio',t]
#            coefkwargs={}
    #           coefkwargs['wave']=wave
    #           coefkwargs['refalt']=refalt
#            coefkwargs['method']=method
    #           coefkwargs['energy']=mpl.header['energy'].ix[t]
    #           coefkwargs['calrange']=calrange
#            if np.isnan(tempprof).any() or np.isnan(templrat).any():
#                print 'trouble at {0}'.format(t)
            if method=='klett2':
                tempbeta,tempsigma=itools.klett2(tempprof,templrat,r_m=refalt,wave=wave)
            tempbackscatter.ix[t]=tempbeta
            tempextinction.ix[t]=tempsigma
    LNCin.backscatter=tempbackscatter
    LNCin.extinction=tempextinction
        
    return LNCin

def Tmatcher(lrat,profin,Tin,method,layerrange,wave,energy,refalt,calrange):    
    alts=profin.index
    if refalt is None:
        for a in alts[::-1]:
            if not np.isnan(profin.ix[a]):
                refalt=a
                break
    
    coefkwargs={}
    coefkwargs['wave']=wave
    coefkwargs['refalt']=refalt
    coefkwargs['method']=method
    coefkwargs['calrange']=calrange
    coefkwargs['energy']=energy
    lratprof=pan.Series(data=lrat,index=profin.index)
    backprof,extprof=coefprofs(profin,lratprof,**coefkwargs)    
    Tcalc=Tfromextinction(extprof.ix[layerrange[0]:layerrange[-1]])
    return Tin-Tcalc

def lratmatcher(backscatter_t,lrat_p,wave=532.0):
    #initial lidar ratios are only 
    lrat_m=8.0*np.pi/3.0
    
    alts=backscatter_t.index.values
    backscatter_m=itools.molprof(alts,wave=wave)['beta_R']
    
    backscatter_p=backscatter_t-backscatter_m
    lrat_total=(lrat_p*backscatter_p+lrat_m*backscatter_m)/(backscatter_t)
    
    return lrat_total

    
def simplelayercorrection(layerprof,backprof,extprof,lrat0,**kwargs):
    #profin must be bounded 
    molbelow=kwargs.get('molbelow',None)
    molabove=kwargs.get('molabove',None)
    method=kwargs.get('method','klett2')
    lratrange=kwargs.get('lratrange',[10,100])
    tolerance=kwargs.get('tolerance',0.01)
    numreps=kwargs.get('numreps',100)
    wave=kwargs.get('wave',1064.0)
    average=kwargs.get('average',True)
    numvals=kwargs.get('numvals',10)
    partrat=kwargs.get('partrat',True)
    energy=kwargs.get('energy',1.0)
    #calculate optical depth from profile or boundaries
    if len(molbelow)>0 and len(molabove)>0:
        Tin=Tfrombounds(molbelow,molabove)
    else:
        Tin=Tfromprofile(layerprof,wave=wave,average=average,numvals=numvals)    
    #calculate optical depth from lidar ratio
    layerrange=[layerprof.index[0],layerprof.index[-1]]    
    if method=='klett2':
        refalt=[]
        calrange=[]
        try:
            lratcor=opt.newton(Tmatcher,lrat0,args=(layerprof,Tin,method,layerrange,wave,energy,refalt,calrange),maxiter=numreps,tol=tolerance)
            print "{0} successfully corrected to {1} after < {2} iterations".format(lrat0,lratcor,numreps)
        except RuntimeError:
            lratcor=lrat0
            print "Warning! Failed to converge after {0} iterations".format(numreps)
        coefkwargs={}
        coefkwargs['refalt']=refalt
        coefkwargs['method']=method
        coefkwargs['wave']=wave
        coefkwargs['calrange']=calrange
        coefkwargs['energy']=energy
        if lratrange[0]<=lratcor<=lratrange[1]:
            lratprof=pan.Series(data=lratcor,index=layerprof.index)
            backprof,extprof=coefprofs(layerprof,lratprof,**coefkwargs)        
        else:
            lratprof=pan.Series(data=lrat0,index=layerprof.index)
            backprof,extprof=coefprofs(layerprof,lratprof,**coefkwargs) 
        if partrat:
                lratout=lratpart(backprof,lratcor)
                lratprof=pan.Series(data=lratout,index=layerprof.index) 
    else:
        tempprof=layerprof.append(molabove)
        refalt=tempprof.index[-1]
        if len(molabove)>numvals:
            calrange=[tempprof.index[-numvals],tempprof.index[-1]]
        else:
            calrange=[molabove.index[0],molabove.index[-1]]
        try:
            lratcor=opt.newton(Tmatcher,lrat0,args=(tempprof,Tin,method,layerrange,wave,energy,refalt,calrange),maxiter=numreps,tol=tolerance)
            print "{0} successfully corrected to {1} after < {2} iterations".format(lrat0,lratcor,numreps)
        except RuntimeError:
            lratcor=lrat0
            print "Warning! Failed to converge after {0} iterations".format(numreps)
        coefkwargs={}
        coefkwargs['refalt']=refalt
        coefkwargs['method']=method
        coefkwargs['wave']=wave
        coefkwargs['calrange']=calrange
        coefkwargs['energy']=energy
        if lratrange[0]<=lratcor<=lratrange[1]:
            lratout=lratcor
        else:
            lratout=lrat0
        molrat=0.0
        lratprof_ref=pan.Series(data=molrat,index=molabove.index)
        lratprof=pan.Series(data=lratout,index=layerprof.index)
        templratprof=lratprof.append(lratprof_ref)
        tempbackprof,tempextprof=coefprofs(tempprof,templratprof,**coefkwargs)
        backprof=tempbackprof.ix[layerprof.index[0]:layerprof.index[-1]]
        extprof=tempextprof.ix[layerprof.index[0]:layerprof.index[-1]]
        if partrat:
            lratprof=lratpart(backprof,lratout)
    
    return backprof,extprof,lratprof

def sceneoptimize(LNCin,**kwargs):
    wave=kwargs.get('wave',1064.0)
    refalt=kwargs.get('refalt',None)
    method=kwargs.get('method','klett2')
    calrange=kwargs.get('calrange',None)
    lrat=kwargs.get('lrat',None)
    inplace=kwargs.get('inplace',True)
    mode=kwargs.get('mode','copol')
        
    if inplace:
        LNCout=LNCin
    else:
        LNCout=deepcopy(LNCin)
    
    #start by performing initial backscatter and extinction corrections using initial lrat values
    coefkwargs={}
    coefkwargs['refalt']=refalt
    coefkwargs['method']=method
    coefkwargs['wave']=wave
    coefkwargs['calrange']=calrange
    coefkwargs['lrat']=lrat
    coefkwargs['mode']=mode
    
    if LNCout.backscatter is None or LNCout.extinction is None:
        LNCout=basiccorrection(LNCout,**coefkwargs)

    BR=LNCout.BR
    backscatter=LNCout.backscatter
    extinction=LNCout.extinction
    panelout=LNCout.scenepanel
        #then go through profiles one by one and perform corrections for each layer


    times=BR.index
    alts=BR.columns
    dflrat=panelout.loc['Lidar_Ratio']
    for t in times:
        BRprof=BR.loc[t]
        lratprof=panelout.loc['Lidar_Ratio',t]
        baseprof=panelout.loc['Base',t]
        topprof=panelout.loc['Top',t]
        typeprof=panelout.loc['Type',t]
        backprof=backscatter.loc[t]
        extprof=extinction.loc[t]
        a=alts[0]
        while True:
            layerbase=baseprof[a]
            layertop=topprof[a]
            lrat=lratprof[a]
            ltype=typeprof[a]
            lratbelow=lratprof[lratprof.index<layerbase]
            lratabove=lratprof[lratprof.index>layertop]
            typebelow=typeprof[typeprof.index<layerbase]
            typeabove=typeprof[typeprof.index>layertop]
            if len(lratabove)==0:
                break
            if len(lratbelow)==0:
                a=lratabove.index[0]
                continue
            if ltype=='PBL' or ltype=='Clear Air' or ltype=='Unidentified Aerosol':
                a=lratabove.index[0]
            else:                
                if typebelow.iloc[-1]=='Clear Air' and typeabove.iloc[0]=='Clear Air':
                    tempBR=BRprof[(BRprof.index>=layerbase)&(BRprof.index<=layertop)]
                    tempback=backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]
                    tempext=extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]
                    tempprof=[]
                    z=[]
                    for tempalt in lratbelow.index[::-1]:
                        if typebelow.ix[tempalt]=='Clear Air':
                            z.append(tempalt)
                            tempprof.append(BRprof.ix[tempalt])
                        else:
                            break
                    molbelow=pan.Series(tempprof[::-1],index=z[::-1])
                    tempprof=[]
                    z=[]
                    for tempalt in lratabove.index:
                        if typeabove.ix[tempalt]=='Clear Air':
                            z.append(tempalt)
                            tempprof.append(BRprof.ix[tempalt])
                        else:
                            break
                    molabove=pan.Series(tempprof,index=z)
                    corkwargs={}
                    corkwargs['energy']=LNCin.header['energy'].ix[t]
                    corkwargs['method']=method
                    corkwargs['molbelow']=molbelow
                    corkwargs['molabove']=molabove
                    
                    print "Type is {0} mean signal is {1}".format(ltype,np.mean(tempBR))
                    backcor,extcor,lratcor=simplelayercorrection(tempBR,tempback,tempext,lrat,**corkwargs)
                    backprof[(backprof.index>=layerbase)&(backprof.index<=layertop)]=backcor
                    extprof[(extprof.index>=layerbase)&(extprof.index<=layertop)]=extcor
                    lratprof[(lratprof.index>=layerbase)&(lratprof.index<=layertop)]=lratcor
            
                a=lratabove.index[0]
            backscatter.loc[t]=backprof
            extinction.loc[t]=extprof
            dflrat.loc[t]=lratprof
        panelout.loc['Lidar_Ratio']=dflrat
    
    LNCout.backscatter=backscatter
    LNCout.extinction=extinction
    LNCout.scenepanel=panelout
               
    return mplout

def AODcalc(mplin):
    total_extinction=mplin.extinction[0]
    typedat=mplin.scenepanel[0]['Type']
    PBLalt=(typedat=='PBL').idxmin(axis=1)
    
    #only include regions of aerosol or PBL
    
    filtered_extinction=total_extinction[(typedat=='Aerosol') | (typedat=='PBL')]
    filtered_extinction.fillna(0.0,inplace=True)
    AOD=pan.Series(index=filtered_extinction.index)
    for t in filtered_extinction.index:
        Ttemp=1.0   #initialize transmissivity to 1.0
        extprof=filtered_extinction.ix[t]  #profile of extinction coefficients in km^-1
        minalt=extprof.index[0]
        #estimate AOD below minalt assuming extinction is same as mean for PBL
        PBLext=np.mean(extprof.ix[:PBLalt.ix[t]])
        Tbelow = np.exp(-minalt*PBLext)        
        Ttemp*=Tbelow
        for lowalt,highalt in zip(extprof.index[:-1],extprof.index[1:]):
            altstep=highalt-lowalt                
            Tstep=np.exp(-altstep*(extprof.ix[lowalt]+extprof.ix[highalt])/2.0) #use trapezoidal integration to calculate step change in transmissivity
            Ttemp*=Tstep 
        AOD.ix[t]=1.0-Ttemp  #AOD is 1.0-T if T is transmissivity
    return AOD
           
        
def quickplot(df,**kwargs):
    ar=kwargs.get('ar',2.0)
    vrange=kwargs.get('vrange',[0,1])
    fsize=kwargs.get('fsize',21)
    maptype=kwargs.get('maptype','customjet')
    orientation=kwargs.get('orientation','Vertical')
    overcolor=kwargs.get('overcolor','w')
    undercolor=kwargs.get('undercolor','k')
    numvals=kwargs.get('numvals',50)
    
    numfigs=len(plt.get_fignums())
    fig=plt.figure(numfigs+1)
    ax=fig.add_subplot(111)    
    data=df.values.T[::-1]
    xdata=df.index
    ydata=df.columns[::-1]
    
    my_cmap=mplot.custom_cmap(maptype=maptype,numvals=numvals,overcolor=overcolor,
                              undercolor=undercolor)  
    im = ax.imshow(data,vmin=vrange[0],vmax=vrange[1],cmap = my_cmap)
    mplot.forceAspect(ax,ar)       
    mplot.altticks(ax, ydata, fsize = fsize, tcolor = 'w')
    
    if orientation=='vertical':
        ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)
    elif orientation=='vorizontal':
        ax.set_ylabel('Horizontal Range [m]', fontsize = fsize+4, labelpad = 15)
    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)        
    ax.axis('tight')
    vstep=(vrange[1]-vrange[0])/5.0
    cbar = fig.colorbar(im, orientation = 'vertical', aspect = 6, extend='both')
    cbar.set_ticks(np.arange(vrange[0],vrange[1]+vstep,vstep))
    cbar.set_ticklabels(np.arange(vrange[0],vrange[1]+vstep,vstep))
    cbar.ax.tick_params(labelsize = fsize)
    cbar.ax.set_ylabel('$[counts*km^{2}/(\mu s*\mu J)$')
    mplot.dateticks(ax, xdata, hours = ['03','06', '09','12', '15','18','21'], 
              fsize = fsize, tcolor = 'w')
    fig.canvas.draw()
    

    
if __name__=='__main__':
    os.chdir('E:\CORALNet\ASCII_Files\Smoke2012\UBC\July')
    savefilepath='E:\CORALNet\ASCII_Files\Smoke2012\UBC\July\Processed'
    plotfilepath='E:\CORALNet\ASCII_Files\Smoke2012\UBC\July\Figures'
    filepath=ltools.get_files('Select HDF5 file',filetype=('.h5','*.h5'))[0]
    filename=os.path.split(filepath)[-1]
    
#    BRfilepath=ltools.get_files('Select BR file',filetype=('.p','*.p'))[0]
#    PRfilepath=ltools.get_files('Select PR file',filetype=('.p','*.p'))[0]
#    MSKfilepath=ltools.get_files('Select MSK file',filetype=('.p','*.p'))[0]

    savefilename='{0}_proc-v3.h5'.format(filename.split('_proc')[0])
    plotfilename='{0}_coefplot.png'.format(filename.split('_proc')[0])
    colorplotfilename='{0}_colorplot.png'.format(filename.split('_proc')[0])
    timestep='600S'
    altrange = np.arange(0.150,15.030,0.030)
    SNRthreshold=1.0
    molthresh=1.0
    layernoisethresh=1.0
    sigma0=0.001
    depolsigma0=0.05
    cloudthresh=(20.0,5.0)
    waterthresh=0.10
    icethresh=0.35
    smokethresh=0.10
    dustthresh=0.20
    maxaeroalt=10.0
    
    PBLCWTrange=np.arange(2,10)
    PBLwidth=2
    doPBL=False
    
    PBLrat=30.0
    molrat=0.0
    moldepol=0.0035
    
    recalc=True
    
    
    print 'Testing LNC process functions'
    
    
    LNCtest = ltools.LNC()
    LNCtest.fromHDF(filepath)
    if recalc:
        LNCtest.alt_resample(altrange)
        LNCtest.calc_all()
        
        layerkwargs={'LNCin':LNCtest,'timestep':timestep,'molthresh':molthresh,
                     'noisethresh':layernoisethresh,'cloudthresh':cloudthresh,
                     'doPBL':doPBL,'PBLCWTrange':PBLCWTrange,'PBLwdith':PBLwidth,
                     'sigma0':sigma0,'depolsigma0':depolsigma0,'waterthresh':waterthresh,
                     'icethresh':icethresh,'smokethresh':smokethresh,'dustthresh':dustthresh,
                     'maxaeroalt':maxaeroalt}
                     
        scenekwargs={'PBLrat':PBLrat,'molrat':molrat,'moldepol':moldepol,'cloudthresh':cloudthresh,
                     'waterthresh':waterthresh,'icethresh':icethresh,'smokethresh':smokethresh,
                     'dustthresh':dustthresh,'maxaeroalt':maxaeroalt}
        layerdict=findalllayers(**layerkwargs)
        LNCtest=scenemaker(layerdict)
        LNCtest=basiccorrection(LNCtest,inplace=True)
        if os.path.isdir(savefilepath):
            os.chdir(savefilepath)
        else:
            os.mkdir(savefilepath)
            os.chdir(savefilepath)
            
        LNCtest.save_to_HDF(savefilename)
    starttime = None
    endtime = None
    hours = []
    topplot_limits=(0.0,5e-5,1e-5)
    bottomplot_limits=(0.0,2e-3,5e-4)
    
    kwargs = {'saveplot':True,'showplot':True,'verbose':True,
                'savefilepath':plotfilepath,'savefilename':plotfilename,
                'hours':hours,'bottomplot_limits':bottomplot_limits,
                'topplot_limits':topplot_limits,'SNRmask':False,
                'colormap':'customjet','orientation':'vertical','toptype':'backscatter',
                'bottomtype':'extinction','logplot':False}
    
    lplot.doubleplot(LNCtest,**kwargs)    
    #   
    colorkwargs = {'saveplot':True,'showplot':True,'verbose':True,
                    'savefilepath':plotfilepath,'savefilename':colorplotfilename,
                    'hours':hours,'SNRmask':False,}
    lplot.colormask_plot(LNCtest,**colorkwargs)


