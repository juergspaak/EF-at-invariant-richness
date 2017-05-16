"""
@author: J.W. Spaak
This files plots Fig. S2
"""
import matplotlib.pyplot as plt
import numpy as np

import help_functions as hef
import community_construction_coex as coex

sqrt = np.sqrt(3)
#get the parameters
t_values = {'e>0,t_e': [par[1] for par in coex.para['e>0']],
        'e<0,t_e': [par[1] for par in coex.para['e<0']],
        'e>0,t_mu': [par[2] for par in coex.para['e>0']],
        'e<0,t_mu': [par[2] for par in coex.para['e<0']],
        'uniform': np.linspace(-1/sqrt, 1/sqrt,len(coex.para['e>0']))}
n = {'e>0,': [par[-1] for par in coex.para['e>0']],
     'e<0,': [par[-1] for par in coex.para['e<0']],
     'uniform': sorted(np.random.randint(5,21,len(coex.para['e>0'])))}

#plotting
fig, (ax1,ax2) = plt.subplots(2,1,figsize =(9,10))
hef.percentiles(t_values, ['t_e', 't_mu'],["red","green"], 
                         plot = (fig,ax1))
ax1.plot(np.linspace(0,100,len(coex.para['e>0'])),t_values['uniform'])

simArtist = plt.Line2D((0,1),(0,0), color='k', marker='v', linestyle='')
anyArtist = plt.Line2D((0,2),(0,0), color='k',marker='^', linestyle='')
handles = []
for col in ['red','green','blue','magenta']:
    handles.append( plt.Line2D((0,2),(0,0), color=col))
ax1.legend(loc = "upper left")
ax1.legend(handles[:3]+[simArtist,anyArtist],
           ['t_e', 't_mu', 'uniform','e<0', 'e>0']
           ,loc = "upper left", numpoints = 1)
hef.percentiles(n, [''],['magenta'], 
                         plot = (fig,ax2))
ax2.plot(np.linspace(0,100,len(coex.para['e>0'])),n['uniform'], )
ax2.set_ylim(ymin = 4, ymax = 21)
ax2.legend(handles[-2:]+[simArtist,anyArtist],
           ['uniform', 'n','e<0', 'e>0']
          ,loc = "upper left", numpoints = 1)
ax1.set_ylabel(r'$t$', fontsize=16)
ax2.set_ylabel(r'$n$', fontsize=16)

#save figure
fig.savefig("Figure S2, Distributions.pdf")