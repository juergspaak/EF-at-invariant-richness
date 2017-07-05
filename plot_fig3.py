"""
@author: J.W. Spaak
This programm plots Fig. 3
"""
import percentiles
import community_construction_repl as repl
import community_construction_coex as coex
  
keys = ['0.95', '0.75', '0.50','0.25', '0.05']

#compute DeltaEF/EF for the different communities
EF_data = {key: repl.delta_EF_lin(*repl.para[key]) for key in repl.para.keys()}
       
# compare to linear functioning        
EF_data["e<0,ref"] = coex.delta_EF_lin(*coex.para["e<0"])
EF_data["e>0,ref"] = coex.delta_EF_lin(*coex.para["e>0"])      
keys = ['ref']+keys
# plot results
ticks = [-80,-60,-40,-20,0,20,40,60,80,100]
labels = ["1.00", "0.95", '0.75',"0.50", 
                      '0.25',"0.05"]

fig, ax,ind = percentiles.bars(EF_data, keys)
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_ylim(-80,80)

ax.set_xticks(ind+0.15)
ax.set_xticklabels(labels)

ax.set_xlabel("proportion p of species present at both sites", fontsize = 16)
fig.savefig("Figure 3, changing community structure, 7barplots.pdf")

