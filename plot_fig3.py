"""
@author: J.W. Spaak
This programm plots Fig. 3
"""
from plot_functions import bars
import community_construction_repl as repl
import community_construction_coex as coex
 
#compute DeltaEF/EF for the different communities
EF_data = {key: repl.delta_EF_lin(*repl.para[key]) for key in repl.para.keys()}
       
# compare to linear functioning        
EF_data["e<0,1.00"] = coex.delta_EF_lin(*coex.para["e<0"])
EF_data["e>0,1.00"] = coex.delta_EF_lin(*coex.para["e>0"])      

# plot results
keys = ['1.00', '0.95', '0.75', '0.50','0.25', '0.05']
fig, ax,ind = bars(EF_data, keys)

# adjust x and y axis
ax.set_ylim(-80,80)
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys)
ax.set_xlabel("Proportion p of species present at both sites", fontsize = 16)

# save the figure
fig.savefig("Figure 3, changing community structure.pdf")

