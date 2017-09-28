"""
@author: J.W. Spaak
This program contains fucntions for plotting
"""
import matplotlib.pyplot as plt
import numpy as np

def bars(datas,keys, fig = None, ax = None, col = ['green', 'red'],
         leg = 'upper right'):
    """ plots boxplots for all data sets in datas
    
    datas is a dictionary with the keys :"e<0,"+key[i] and "e>0,"+key[i],
        the 0.05, 0.25, 0.50, 0.75, 0.95 percentiles of these will be plotted
    fig, ax: figure and axis of the figure if boxes should be plotted
        into a preexisting figure
    keys are the keys of the datas
    color = [color[0], color[1]]
    leg is the location of the legend
    
    returns fig,ax the figure and the axis of the figure
            ind index on x value of the boxes"""
    per_pos = {} # dictionary containing the percentiles of the EF, for e>0
    per_neg = {} # dictionary containing the percentiles of the EF, for e<0
    
    # compute the percentiles for all EF_data
    for per in [5,25,50,75,95]:
        per_pos[str(per)] = []
        per_neg[str(per)] = []
        for key in keys:
            # append the percentiles
            per_pos[str(per)].append(np.percentile(datas['e<0,'+key], per))
            per_neg[str(per)].append(np.percentile(datas['e>0,'+key], per))
        #convert to array
        per_pos[str(per)] = np.array(per_pos[str(per)])
        per_neg[str(per)] = np.array(per_neg[str(per)])
    
    #width of boxes
    width = 0.4
    ind = np.arange(len(keys))*0.7 #index of boxes
    if fig == None: # if plot needs to be continues
        fig, ax = plt.subplots(figsize = (9,7))
        zero_line = True
    else:
        ind += 0.05 #to avoid box overlaping
        zero_line = False
    
    # error bars (5 and 95 percentile)
    yerr = [per_pos['75']-per_pos['5'],per_pos['95']-per_pos['75']]
    #plot the boxes
    pos_bar = ax.bar(ind, per_pos['75']-per_pos['25'], width,  yerr = yerr, 
           bottom = per_pos['25'],edgecolor = col[0],linewidth= 1.5,
           fill =False, ecolor = col[0])
    
    ind += width/4 #shift the boxes a little bit
    # error bars (5 and 95 percentile)
    yerr = [per_neg['75']-per_neg['5'],per_neg['95']-per_neg['75']]
    #plot the boxes
    neg_bar = ax.bar(ind, per_neg['75']-per_neg['25'], width, yerr = yerr,
           bottom = per_neg['25'], edgecolor = col[1], linewidth= 1.5, 
           fill =False, ecolor = col[1])
    
    # plot the medians
    plt.scatter(ind+width/4, per_pos['50'], 40,col[0])
    plt.scatter(ind+width/2, per_neg['50'], 40,col[1])
    
    # adjust axis
    ax.set_yticks([-80,-60,-40,-20,0,20,40,60,80,100])
    ax.set_ylim(-80,100)
    ax.set_ylabel(r'$100\cdot\Delta EF/EF_u$', fontsize=16)
    
    # add legend
    ax.legend([pos_bar[0],neg_bar[0]],['e<0', 'e>0'], loc = leg)
    if zero_line:# add zero line
        ax.plot([-0.2, ind[-1]+0.8], [0,0], '--', color = "0.5")
    return fig, ax, ind

 
def percentiles(datas,keys,color, y_min=None, y_max = None,labels = None,
                       fsize = 16,ticks = False, ls = None, plot = None):
    """ plots the percentile curves for all data sets in datas
    
    datas is a dictionary with the keys :"e<0,"+key[i] and "e>0,"+key[i],
        the percentile curves these will be plottet
    keys are the keys of the datas
    color = [color[0], color[1],...]
    y_min and y_max are the axis settings
    labels are shown in legend, by default labels = keys
    fsize is the font size
    ticks are the ticks on the y axis
    ls is the line style
    plot =(fig,ax) can be used to continue plotting on a given plot
    
    returns fig,ax, the figure and the axis of the figure"""
    datas_pos  = {key: datas["e>0,"+key] for key in keys} #split into two cases
    datas_neg  = {key: datas["e<0,"+key] for key in keys}

    if plot is None: #start new plot
        fig, ax = plt.subplots(figsize = (9,7))
    else: #continues plot on previous plot
        fig,ax = plot
    y_Min, y_Max = 0,0
    acc = 1000 #number of point plotted
    if labels is None: #only used for fast/exploratory work
        labels = keys
    if ls is None: #only used for fast/exploratory work
        ls = len(keys)*['-']
    
    for iteration, key in list(enumerate(keys)):
        data_pos_pre = sorted(datas_pos[key].copy()) #sort for percentile curves
        data_neg_pre = sorted(datas_neg[key].copy())
        data_pos = data_pos_pre[::int(len(data_pos_pre)/acc)] #sample 1000 points
        percent_pos = np.linspace(0,100,len(data_pos),endpoint = False)
        data_neg = data_neg_pre[::int(len(data_neg_pre)/acc)]
        percent_neg = np.linspace(0,100,len(data_neg),endpoint = False)
        # if len(data_pos_pre)%1000 !=1, then we miss the maximum
        data_pos = np.append(data_pos,data_pos_pre[-1]) #append maximum
        percent_pos = np.append(percent_pos,100)
        data_neg = np.append(data_neg,data_neg_pre[-1]) #append maximum
        percent_neg = np.append(percent_neg,100)
        #plot curves
        ax.plot(percent_pos,data_pos, color = color[iteration], 
                linestyle = ls[iteration],  marker = '^', 
                markevery = 50, label = labels[iteration])
        ax.plot(percent_neg,data_neg, color = color[iteration], 
                linestyle = ls[iteration],  marker = 'v', 
                markevery = 50, label = labels[iteration])
        y_Min = min(y_Min, data_pos[0], data_neg[0]) #min and max for plot
        y_Max = max(y_Max, data_pos[-1],data_neg[-1])
    #axis labels
    fig.gca().set_ylabel(r'$100\cdot\Delta EF/EF_u$', fontsize=fsize)
    fig.gca().set_xlabel("percentile", fontsize=fsize)
    #min and max predefined by user?
    if y_min is None:
        y_min = y_Min
    if y_max is None:
        y_max = y_Max    
    ax.axis([0,100,y_min, y_max]) #set min and max of plot
    ax.tick_params(axis='both', which='major', labelsize=fsize-2)
    if ticks:
        plt.yticks(ticks)
    #legend settings
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='v', linestyle='')
    anyArtist = plt.Line2D((0,2),(0,0), color='k',marker='^', linestyle='')
    handles = []

    for i in range(len(keys)):
        handles.append( plt.Line2D((0,2),(0,0), color=color[i], ls = ls[i]))
    ax.legend(loc = "upper left")
    ax.legend(handles+[simArtist,anyArtist],
              labels+['e<0', 'e>0']
               ,loc = "upper left", numpoints = 1)
    return fig,ax  

