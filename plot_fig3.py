"""
@author: J.W. Spaak
This programm plots Fig. 3
"""
import percentiles
import community_construction_repl as repl
import community_construction_coex as coex
  
keys = ['0.95', '0.50', '0.05','rand']
#compute DeltaEF/EF for the different communities
EF_data = {}
for k in keys:
    EF_data["e<0,"+k] = repl.delta_EF_lin(*repl.para["e<0,"+k])
    EF_data["e>0,"+k] = repl.delta_EF_lin(*repl.para["e>0,"+k])
        
# compare to linear functioning        
EF_data["e<0,ref"] = coex.delta_EF_lin(*coex.para["e<0"])
EF_data["e>0,ref"] = coex.delta_EF_lin(*coex.para["e>0"])      
keys = ['ref']+keys
# plot results
ticks = [-80,-60,-40,-20,0,20,40,60,80,100]
labels = ["p = 1.00", "p = 0.95", "p = 0.50", "p = 0.05","p is random"]
blues = 0*np.ones([5,3])
blues[:,0] = np.linspace(0,0.5,5)[::-1]

blues[:,-1] = np.linspace(.5,1,5)


fig, ax,ind = percentiles.bars(EF_data, keys, color)
ax.set_xlim([-0.2,3.5])

ax.set_xticks(ind+0.15)
ax.set_xticklabels(labels)

ax.set_xlabel("Percent of species present at both sites", fontsize = 16)
fig.savefig("Figure 3, changing community structure, barplots.pdf")

#plottig percentile curves
fig,ax = percentiles.percentiles(EF_data, keys,y_min = -80, y_max = 100, 
        ticks = ticks,color = color, labels = labels)


#save figure
fig.savefig("Figure 3, changing community structure.pdf")