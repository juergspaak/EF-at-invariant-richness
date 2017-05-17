"""
@author: J.W. Spaak
This program was used to make Fig.5
"""
import numpy as np
import matplotlib.pyplot as plt

import help_functions as hef
import community_construction_coex as coex

#different values for H that are used:
keys = ['100','2', '1','0.5','0.1']

fig, (ax1,ax2) = plt.subplots(1,2,figsize =(18,7))
# asymptotic functions used
H = lambda N, aveH: aveH*N/(aveH+N)
# plot these different functions
color = ["red", "green", "blue", "black", "magenta","orange"]
H_functions = {key: H(np.linspace(0,1,1000),float(key)) for key in keys[:5]}
for (i,key) in enumerate(keys):
    ax1.plot(np.linspace(0,1,1000),H_functions[key], label = "H = "+key
             ,color = color[i])

# x-axis: densitiy of species relative to equilibrium density
ax1.set_xlabel(r'$N/N^*$', fontsize=16)
# y-axis: total fun contributed by this species (not per capita contribution)
ax1.set_ylabel(r'$f(N)$', fontsize=16)
ax1.legend(loc = "upper left", fontsize = 16)
ax1.set_title("A. Examples of asymptotic functions", fontsize = 16)

keys = ['100','2', '1','0.5','0.1']
# contains all EF datas
EF_data = {}
#compute asymptotic functioning for different averages of H
for key in keys:
    EF_data["e<0,"+key] = \
        hef.comp_EF(coex.delta_EF_asym,coex.para["e<0"],max_ave_H=float(key))
    print("Progress report: e<0, H="+key+" is done" )
    EF_data["e>0,"+key] = \
        hef.comp_EF(coex.delta_EF_asym,coex.para["e>0"],max_ave_H=float(key))
    print("Progress report: e>0, H="+key+" is done" )
    
#compare to linear functioning        
EF_data["e<0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,coex.para["e<0"])
EF_data["e>0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,coex.para["e>0"])
keys.append("ref")
# plot results
ticks = [-60,-40,-20,0,20,40,60,80,100]

labels = [r'max$(\overline{H})=100$',r'max$(\overline{H})=2$',
          r'max$(\overline{H})=1$',r'max$(\overline{H})=0.5$',
          r'max$(\overline{H})=0.1$',"Ref. Linear functioning"]
fig,ax2 = hef.percentiles(EF_data, keys, y_max = 100, ticks = ticks,
       labels = labels,color = color, plot = (fig,ax2))
ax2.set_title('B. EF with asymptotic functioning', fontsize = 16)

#save figure
fig.savefig("Figure 4, Asymptotic functioning.pdf")