"""
@author: J.W. Spaak
This program contains fucntions for plotting and computing several communities
"""
import matplotlib.pyplot as plt
import numpy as np
import community_construction_coex as coex


def com_con(rand_par, num_com, **kwargs):
    com_paras = []
    for i in range(num_com):
        com_paras.append(rand_par(**kwargs))
    return com_paras
    
def comp_EF(delta_EF, com_para,**kwargs):
    EF_data = np.zeros(len(com_para))
    for (i,para) in list(enumerate(com_para)):
        EF_data[i] = delta_EF(*com_para[i],**kwargs)
    return EF_data
    
def percentiles(datas,keys, y_min=None, y_max = None,labels = None,
                     color = None,  fsize = 16,ticks = False):
    """ plots the percentile curves for all data sets in data
    
    datas = [data1, data2...]: is a list with lists, each list will be plotted
    or datas can be a dictionary, then labels are assumed (unless otherwise specified)
    to be the name tags in the dictionary
    labels are the labels of the datas in the legend, must be ordered similarly
    y_label is the labe on the y_label (on y_axis)
    start, end are the values on x axis
    acc is the accuracy of the plot
    y_min, y_max are the boundaries of the y_axis
    save = True saves the plot with y_label name, alternatively with save, then save must be a string
    plots contains instructions to plot each line
    ordered_keys can be used to have different order in the legend (Default is dictionary order)"""
    datas_pos  = {key: datas["e>0,"+key] for key in keys}
    datas_neg  = {key: datas["e<0,"+key] for key in keys}
    fig, ax = plt.subplots(figsize = (9,7))
    y_Min, y_Max = 0,0
    acc = 1000
    if labels == None:
        labels = keys
    
    for iteration, key in list(enumerate(keys)):
        data_pos = sorted(datas_pos[key].copy())
        data_neg = sorted(datas_neg[key].copy())
        data_pos = data_pos[::int(len(data_pos)/acc)]
        percent_pos = np.linspace(0,100,len(data_pos),endpoint = False)
        data_neg = data_neg[::int(len(data_neg)/acc)]
        
        percent_neg = np.linspace(0,100,len(data_neg),endpoint = False)
        plt.plot(percent_pos,data_pos, color = color[iteration], 
                marker = '^', markevery = 50, label = labels[iteration])
        plt.plot(percent_neg,data_neg, color = color[iteration], 
                marker = 'v', markevery = 50,label = labels[iteration])
        y_Min = min(y_Min, data_pos[0], data_neg[0])
        y_Max = max(y_Max, data_pos[-1],data_neg[-1])
    
    fig.gca().set_ylabel(r'$100\cdot\Delta EF/EF$', fontsize=fsize)
    fig.gca().set_xlabel("percentile", fontsize=fsize)
    if y_min is None:
        y_min = y_Min
    if y_max is None:
        y_max = y_Max    
    ax.axis([0,100,y_min, y_max])
    plt.tick_params(axis='both', which='major', labelsize=fsize-2)
    if ticks:
        plt.yticks(ticks)

    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='v', linestyle='')
    anyArtist = plt.Line2D((0,2),(0,0), color='k',marker='^', linestyle='')
    handles = [1,2,3,4,5]
    display = range(5)
    for i in range(5):
        handles[i] = plt.Line2D((0,2),(0,0), color=color[i])
    ax.legend(loc = "upper left")
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[simArtist,anyArtist],
              labels+['e<0', 'e>0']
               ,loc = "upper left", numpoints = 1)
    return fig,ax
    

    

