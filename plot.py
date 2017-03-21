import numpy as np
import matplotlib.pyplot as plt
from timeit import default_timer as timer
"""computes and plots the number of species for a maximal EF,
depending on rho, the covariance and alpha, the interaction coefficient"""



def plotter(x_min,x_max,function,accuracy=1000, ls = '-'):
    """ plots the function between x_min, x_max"""
    x_values=np.linspace(x_min,x_max,accuracy)
    
    y_values=[function(x) for x in x_values]
    plt.plot(x_values,y_values,ls)
    plt.grid()
    
def plit(xdata,ydata, c=[None],axis=None,x_label=None, y_label=None,
         xscale='linear', yscale='linear', grid = True):
    """plots the parameters a,b and their c=color value"""
    if len(c) != len(xdata):
        c=np.ones(len(xdata))
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.set(xscale=xscale, yscale=yscale)
    plt.scatter(xdata,ydata,c=c,lw=0,s=1)
    
    if axis: ax.axis(axis)
    if np.var(c)>1e-5: plt.colorbar()
    if grid: plt.grid()
    
    if x_label:  fig.gca().set_xlabel(x_label, fontsize=14)   
    if y_label: fig.gca().set_ylabel(y_label, fontsize=14)
    
    
def plot_all(data, colors=[None], names=None, axis=None, logs=None):
    iterations = len(data)
    """calls plit for all combination of data"""
    if names is None:
        names = range(iterations)
        
    if logs is None:
        logs = iterations*['linear']
        
    if axis is None:  
        for i in range(iterations):
            for j in range(iterations-i-1):
                plit(data[i],data[i+j+1], colors, axis, names[i], names[i+j+1],
                     logs[i], logs[i+j+1])
    else:
        for i in range(iterations):
            for j in range(iterations-i-1):
                plit(data[i],data[i+j+1], colors, axis[i]+axis[i+j+1],
                     names[i], names[i+j+1],logs[i], logs[i+j+1])

def plot_percentile(real_data, ylabel=None, a=0,b=100, acc=1000, save = None):
    """plots the percentiles of real_data between a and b"""
    datas = real_data.copy()
    datas.sort()
    percent=np.zeros(acc)
    t=np.linspace(a,b,acc)
    for i in range(acc):
        percentile = len(datas)/100*(a+(b-a)*i/acc)
        percent[i]=datas[int(percentile)]
    fig, ax = plt.subplots(figsize = (9,7))
    plt.plot(t,percent)
    fig.gca().set_xlabel("percentiles", fontsize=14)
    fig.gca().set_ylabel(ylabel, fontsize = 14)
    plt.grid()
    if save:
        if save is True:
            save_fig("percentile "+ylabel, fig)
        else:
            save_fig("percentile "+save, fig)
    return fig

    
def plot_percentiles(datas,labels = None, y_label=None, start=0,end=100, acc=1000,
                     y_min=None, y_max = None, save = None, grid = True,
                     color = None, ls = None, ordered_keys = None, fsize = 16):
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
    if ordered_keys is None:
        try:
            keys = list(sorted(datas.keys()))
            
        except AttributeError:
            keys = range(len(datas))
    else:
        keys = ordered_keys
       
    if labels is None:
            labels = keys
    if color is None:
        color = len(keys)*[None]
    if ls is None:
        ls = len(keys)*[None]
        
    fig, ax = plt.subplots(figsize = (9,7))
    y_Min, y_Max = 0,0
    
    for iteration, key in list(enumerate(keys)):
        data_i = datas[key].copy()
        data = sorted(data_i)
        
        percent = np.zeros(acc)
        t=np.linspace(start,end,acc)
        for j in range(acc):
            percentile = len(data)/100*(start+(end-start)*j/acc)
            percent[j]=data[int(percentile)]
        plt.plot(t,percent, color = color[iteration], 
                 linestyle = ls[iteration], label=labels[iteration])
        y_Min = min(y_Min, data[0])
        y_Max = max(y_Max, data[-1])
     
    if y_label:
        fig.gca().set_ylabel(y_label, fontsize=fsize)
    fig.gca().set_xlabel("percentile", fontsize=fsize)
    if y_min is None:
        y_min = y_Min
    if y_max is None:
        y_max = y_Max
    ax.axis([start,end,y_min, y_max])
    plt.tick_params(axis='both', which='major', labelsize=fsize-2)
    if grid: plt.grid()
    plt.legend(labels,loc="upper left", fontsize = fsize)
    if save:
        if save is True:
            save_fig("percentiles "+y_label, fig)
        else:
            save_fig("percentiles "+save, fig)
    return fig

def save_fig(save, fig, typ = "pdf"):
    fig.savefig("Figure, " + save+'.'+typ, format = typ, dpi=fig.dpi) #saves with same resolution
    
def time_it(fun, *args, iterations = 10):
    if args == ():
        start = timer()
        for i in range(iterations):
            fun()
        end = timer()
    else:
        start = timer()
        for i in range(iterations):
            fun(*args)
        end = timer()
    print("one run takes about",(end-start)/iterations, "seconds")
    
    
print ("The following functions now exist:"+'\n\n'
       "plotter(plotter(x_min,x_max,function,accuracy=1000):"+'\n'
       "  plots the function between x_min, x_max"+'\n\n'
       
       "plit(xdata,ydata, c=[None],axis=None,x_label=None, y_label=None,"
        "xscale='linear', yscale='linear'):" +'\n'
       """  plots the parameters a,b and their c=color value""" +'\n\n'
        
        "plot_all(data, colors=[None], names=None, axis=None, logs=None):"
        "iterations = len(data)"+'\n'
        """  calls plit for all combination of data""" +'\n\n'
        
        "plot_percentile(real_data, ylabel=None, a=0,b=100, acc=1000):"+'\n'
        """plots the percentiles of real_data between a and b"""+'\n\n'
        
        "plot_percentiles(datas,labels = None, ylabel=None, a=0,b=100, acc=1000,"
        "y_min=None, y_max = None):"+'\n'
        """  plots the percentile curves for all data sets in data""")
            
           
    

