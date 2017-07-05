"""
@author: J.W. Spaak
This programm plots Fig. 4
"""
import numpy as np
import matplotlib.pyplot as plt

import percentiles
import community_construction_coex as coex

#different values for H that are used:



keys = ['100','2', '1','0.5','0.1']
# contains all EF datas
EF_data = {}
#compute asymptotic functioning for different averages of H
for key in keys:
    EF_data["e<0,"+key] = coex.delta_EF_asym(\
                    *coex.para["e<0"],max_ave_H=float(key))
    print("Progress report: e<0, H="+key+" is done" )
    EF_data["e>0,"+key] = coex.delta_EF_asym(\
                    *coex.para["e>0"],max_ave_H=float(key))
    print("Progress report: e>0, H="+key+" is done" )
    
#compare to linear functioning        
EF_data["e<0,ref"] = coex.delta_EF_lin(*coex.para["e<0"])
EF_data["e>0,ref"] = coex.delta_EF_lin(*coex.para["e>0"])
keys = ['ref', '100','2', '1','0.5','0.1']
# plot results

labels = ['Ref',r'max$(\overline{H})=100$',r'max$(\overline{H})=2$',
          r'max$(\overline{H})=1$',r'max$(\overline{H})=0.5$',
          r'max$(\overline{H})=0.1$']
fig, ax, ind = percentiles.bars(EF_data, keys)

ax.set_xlim([-0.2,ind[-1]+0.7])

ax.set_xticks(ind+0.15)
ax.set_xticklabels(labels)
ax.set_ylim(-60, 80)
pos = ax.bar(0,0,width = 0 ,color = 'white',edgecolor='red',linewidth=1.5)
neg = ax.bar(0,0,width = 0 ,color = 'white',edgecolor='green',linewidth=1.5)

ax.legend([neg,pos],['e<0', 'e>0'], loc = 'upper right')
fig.savefig("Newfigure S_, asymptotic fucntioning, 5barplots.pdf")