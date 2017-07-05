"""
@author: J.W. Spaak
This programm plots Fig. S5,S6
"""
import numpy as np
from scipy.stats import linregress

import percentiles
import community_construction_coex as coex
import community_construction_repl as repl

# computes the relative change in average growthrate (mu-mu')/mu
eff_ave = {'e>0,1.00': coex.para['e>0'][0],
        'e<0,1.00': coex.para['e<0'][0]}

#same for species replacing
for key in repl.para.keys():
    mu, e,comp,f,p, alpha = repl.para[key]
    mu_av = p*mu['avb']+(1-p)*mu['avu'] #average growth rate in ref site
    eff_ave[key] = (mu_av-mu['av_change'])/mu_av #change in av growth rate
    


#plot S5
keys = ['1.00','0.95', '0.75', '0.50', '0.25', '0.05']
S5,ax,ind = percentiles.bars(eff_ave, keys)

# adjust axis
ax.set_ylabel(r"$(\overline{\mu}-\overline{\mu'})/\overline{\mu}$",
               fontsize = 16)
ax.set_ylim(-0.5, 0.8)
ax.set_yticks(np.linspace(-0.4,0.8,7))
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys)
ax.set_xlabel("Proportion p of species present at both sites", fontsize = 16)
# save S5
S5.savefig("Figure S5, Adapted e.pdf")

############ For Fig S6
print("Now starting to generate new communities")

adapted_para = {}#contains parameter of communities with adapted e distribution
for key in sorted(eff_ave.keys()):
    if key[4:] == 'coex': #don't need to redo these
        continue
    if key == 'e>0,0.05':
        print("Progress report: e<0 communities have been constructed")
    num = len(eff_ave[key]) #number of communities
    #linearly approximate the change in \mu
    slope, intercept, a,b,c = \
        linregress(np.linspace(0,100, num),sorted(eff_ave[key]))
    ave_max = 100*slope+intercept # maximum for average
    ave_min = intercept #minimum for average
    dist = min(ave_min -(-1), 1-ave_max) #distance to boundaray [-1,1]
    e_min, e_max = ave_min-dist, ave_max+dist #to have symmetric conditions
    # construct communities with adapted e distribution
    adapted_para[key] = coex.rand_par(ave_min = ave_min, ave_max = ave_max,
                           e_min = e_min, e_max = e_max, num=num)
    
print("Progress report: e>0 communities have been constructed")    

# compute the EF for the new computed communities    
EF_data = {key: repl.delta_EF_lin(*repl.para[key]) for key in repl.para.keys()}    
EF_data_adapted = {key: coex.delta_EF_lin(*adapted_para[key]) 
                            for key in adapted_para.keys()}

# plot S6                            
cols = ['#CF0000', '#90FB90','#800000', '#006400']
# plot the delta EF of the invariant community structure
S6, ax, ind = percentiles.bars(EF_data, keys[1:], col =cols[:2])
# add the deltaEF of the adapted e communities
S6, ax, ind = percentiles.bars(EF_data_adapted, keys[1:], 
                                S6, ax, col = cols[2:])

# adjust axis
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys[1:])
ax.set_xlabel("Proportion p of Species present at both Sites", fontsize = 16)
ax.set_ylim([-80,100])

# add legend
legend = {col: ax.bar(0,0,width = 0 ,color = 'white'
                      ,edgecolor=col,linewidth=1.5) for col in cols}
lab = ['e>0, inv. com.', 'e<0, inv. com.', 'e>0, var. com.', 'e<0, var. com.']
ax.legend([legend[col] for col in cols],lab, loc = 'upper right')

# save S6
S6.savefig("Figure S6, Adapted e.pdf")