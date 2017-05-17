"""
@author: J.W. Spaak
This program was used to make Fig.S3
"""
import matplotlib.pyplot as plt

import community_construction_repl as repl
import help_functions as hef


keys = ['0.95', '0.50', 'rand', '0.05']

#compute DeltaEF/EF for the different communities
EF_data = {}
for k in keys:
    #compute with adjustment terms
    EF_data["e<0,"+k] = \
        hef.comp_EF(repl.delta_EF_lin,repl.para["e<0,"+k])
    EF_data["e>0,"+k] = \
        hef.comp_EF(repl.delta_EF_lin,repl.para["e>0,"+k])

    #compute with no adjustment terms
    EF_data["e<0,no adjust"+k] = \
        hef.comp_EF(repl.delta_EF_lin,repl.para["e<0,"+k], adjust = 0)
    EF_data["e>0,no adjust"+k] = \
        hef.comp_EF(repl.delta_EF_lin,repl.para["e>0,"+k], adjust =0)

#plotting
fig, ax =  plt.subplots(2,2,figsize =(12,12), sharex = True, sharey = True)
fig.subplots_adjust(hspace=0.1,wspace = 0.1)
ls = [':', '--']
ticks = [-80,-60,-40,-20,0,20,40,60,80,100]
color = ['green','red', 'cyan', 'purple']
for i,key in list(enumerate(keys)):
    pl_keys = ['no adjust'+key, key]
    col_pl = 2*[color[i]]
    ax_pl = ax[i%2,i//2]
    labels = ["no adjustment terms", "with adjustment terms"]
    hef.percentiles(EF_data, pl_keys,y_min = -80, y_max = 100, 
        ticks = ticks,color = col_pl, ls = ls, labels = labels, 
        plot = [fig, ax_pl])
    ax_pl.set_title('p = ' +key, fontsize = 16)
    
ax[1][1].set_ylabel(' ')
ax[1][0].set_xlabel('percentile', fontsize = 16)
ax[0][0].set_ylabel(r'$100\cdot\Delta EF/EF$', fontsize = 16)
ax[1][0].set_ylabel(r'$100\cdot\Delta EF/EF$', fontsize = 16)
ax[0,1].set_title('p is random', fontsize = 16)

#save figure
fig.savefig("Figure S3, Adjustment terms.pdf")
