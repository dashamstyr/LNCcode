def depol_plot(fig, ax, xdata, ydata, data, fsize = 21):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identically zero because they couldn't be calculated
    cdict = {'red': ((0,1,1),
                     (0.0001, 1, 0),
                     (0.35, 0, 0),
                     (0.66, 1, 1),
                     (0.89,1, 1),
                     (1, 0.5, 0.5)),
         'green': ((0,1,1),
                   (0.0001, 1, 0),
                   (0.125,0, 0),
                   (0.375,1, 1),
                   (0.64,1, 1),
                   (0.91,0,0),
                   (1, 0, 0)),
         'blue': ((0,1,1),
                  (0.0001,1,0.5),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65,0, 0),
                  (1, 0, 0))}
    
    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,1064)
    
    im = ax.imshow(data, vmin=0, vmax=0.5, cmap = my_cmap)
    forceAspect(ax,ar)
        
    altticks(ax, ydata, fsize = fsize)

    
    ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)

    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)
        
    ax.axis('tight')

    return im

def backscatter_plot(fig, ax, xdata, ydata, data, fsize = 21):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors
    import numpy as np
    
    #set colormap to be the same as 'jet' with the addition of white color for
    #depol ratios set to identiacally zero because they couldn't be calculated
    cdict = {'red': ((0,0,0),
                     (0.35, 0, 0),
                     (0.66, 1, 1),
                     (0.89,1, 1),
                     (1, 0.5, 0.5)),
         'green': ((0,0,0),
                   (0.125,0, 0),
                   (0.375,1, 1),
                   (0.64,1, 1),
                   (0.91,0,0),
                   (1, 0, 0)),
         'blue': ((0,0.1,0.1),
                  (0.11, 1, 1),
                  (0.34, 1, 1),
                  (0.65,0, 0),
                  (1, 0, 0))}
    
    my_cmap = colors.LinearSegmentedColormap('my_colormap',cdict,1064)
    
    im = ax.imshow(data, vmin=1, vmax=10, cmap = my_cmap)
    forceAspect(ax,ar)       
    altticks(ax, ydata, fsize = fsize, tcolor = 'w')

    
    ax.set_ylabel('Altitude [m]', fontsize = fsize+4, labelpad = 15)

    for line in ax.yaxis.get_ticklines():
        line.set_markersize(10)
        line.set_markeredgewidth(1)
        
    ax.axis('tight')

    return im
    

def forceAspect(ax,aspect=1):
    im = ax.get_images()
    extent =  im[0].get_extent()
    ax.set_aspect(abs((extent[1]-extent[0])/(extent[3]-extent[2]))/aspect)

def dateticks(ax, axisdat,hours = [], fsize = 21, tcolor = 'k'):
    import matplotlib.pyplot as plt
    from time import strftime
    
    dold = axisdat[0].strftime('%d')
    hold = axisdat[0].strftime('%H')
    tickmarks = []
    ticklabels = []
    n = 0
    l = len(axisdat)

    
    for d in axisdat:
        dtemp = d.strftime('%d')
        if dtemp != dold:
            ticklabels.append(d.strftime('%H\n%b %d'))
            tickmarks.append(n)
        else:
            htemp = d.strftime('%H')
            mtemp = d.strftime('%M')
            if not hours:
                if htemp != hold:
                    ticklabels.append(d.strftime('%H'))
                    tickmarks.append(n)
            else:
                if htemp in hours and htemp != hold:
                    ticklabels.append(d.strftime('%H'))
                    tickmarks.append(n)
            hold = htemp

        dold = dtemp
        n += 1
    
    plt.xticks(tickmarks,ticklabels,fontsize = fsize)

    for line in ax.xaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(2)

def altticks(ax, axisdat, numticks = 5, fsize = 21, tcolor = 'k'):
    import matplotlib.pyplot as plt

    numpoints = len(axisdat)
    step = numpoints//numticks
    tickmarks = range(0,numpoints,step)
    ticklabels = [str(int(t)) for t in axisdat[::step]]

    plt.yticks(tickmarks,ticklabels, fontsize = fsize)

    for line in ax.yaxis.get_ticklines():
        line.set_color(tcolor)
        line.set_markersize(10)
        line.set_markeredgewidth(3)
        

    
if __name__ == '__main__':
    import pandas as pan
    import os
    import LNC_tools as LNC
    import matplotlib.pyplot as plt
    import numpy as np
    import datetime as dt

#    os.chdir('K:\CORALNet\Data\ASCII Files')

    #set date(s) for plot title
##    days = list(set(t.strftime('%B %d, %Y') for t in xdata))
##    days.sort()
##    sdays = list(set(t.strftime('%Y%m%d') for t in xdata))
##    sdays.sort()
##    
##    if len(days) == 1:
##        titledays = days[0]
##        savedays = sdays[0]
##    else:
##        titledays = '%s - %s'%(days[0],days[-1])
##        savedays = '%s-%s'%(sdays[0],sdays[-1])
##
##    title = 'Depolarization Ratio: %s'%titledays
##    savetitle = 'Depolrat_%s'%savedays
    
    #create figure and plot image of depolarization ratios
    fsize = 18 #baseline font size
    ar = 2.0  #aspect ratio
    figheight = 12 #inches

    plt.rc('font', family='serif', size=fsize)

##    cbar = fig.colorbar(im, orientation = 'horizontal', pad = 0.15, aspect = 40)
##    cbar.ax.tick_params(labelsize = font)
##
##    t = ax.set_title(title, fontsize = fsize+10)
##    t.set_y(1.1)
##
##    plt.subplots_adjust(top = 0.86, bottom = 0.01, left = 0.09, right = 0.95)
##
    
    startdate = dt.datetime(2013,4,29,17,0)
    enddate = dt.datetime(2013,4,30,16,0)
    minalt = 500
    maxalt = 10000
    
    fig = plt.figure()

    h_set = ['12']
    
    olddir = os.getcwd()
    
    os.chdir('C:\Users\dashamstyr\Desktop\CORALNet - May2013\ASCII')
    
    filepath = LNC.get_files('Select first file to be plotted', filetype = ('.pickle','*.pickle'))
    
    if filepath[0] == '{':
        filepath = filepath[1:-1]
    
    [path, filename] = os.path.split(filepath)
    
    os.chdir(path)
    df = pan.load(filename)
    df = df.loc[startdate:enddate,:maxalt]
    if minalt != 0:
        df.loc[:,:minalt] = 'nan'
    
    datetime = df.index
    alt = df.columns
    ax = fig.add_subplot(2,1,1)
    im = backscatter_plot(fig, ax, datetime,alt[::-1],df.T[::-1], fsize = fsize)
    cbar = fig.colorbar(im, orientation = 'vertical', aspect = 6)
    cbar.ax.tick_params(labelsize = fsize)
    dateticks(ax, datetime, hours = h_set, fsize = fsize, tcolor = 'w')
    ax.set_xticklabels([])

    filename = LNC.get_files('Select second file to be plotted', filetype = ('.pickle','*.pickle'))
    
    if filename[0] == '{':
        filename = filename[1:-1]
        
    df = pan.load(filename)
    df = df.loc[startdate:enddate,:maxalt]
    if minalt != 0:
        df.loc[:,:minalt] = 'nan'
    datetime = df.index
    alt = df.columns
    ax = fig.add_subplot(2,1,2)
    im = depol_plot(fig, ax, datetime,alt[::-1],df.T[::-1], fsize = fsize)
    cbar = fig.colorbar(im, orientation = 'vertical', aspect = 6)
    cbar.set_ticks(np.arange(0,0.6,0.1))
    cbar.set_ticklabels(np.arange(0,0.6,0.1))
    #set axis ranges and tickmarks based on data ranges
    dateticks(ax, datetime, hours = h_set, fsize = fsize)
    ax.set_xlabel('Time [PDT]',fontsize = fsize+4)
    fig.autofmt_xdate()
        

        ##plt.savefig(savetitle,dpi = 100, edgecolor = 'b', bbox_inches = 'tight')
    fig.set_size_inches(figheight*ar,figheight) 
    fig.tight_layout()
    plt.show() 
    
    os.chdir(olddir)
    
