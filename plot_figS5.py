"""
@author: J.W. Spaak
This programm plots Fig. S4
"""
from plot_functions import bars
import community_construction_repl as repl
  


#compute DeltaEF/EF for the different communities and different cases
EF_data = {key: repl.delta_EF_lin(*repl.para[key]) for key in repl.para.keys()}
# species have different f
EF_data_dif = {key: repl.delta_EF_lin(*repl.para[key], sim_f = False)
                 for key in repl.para.keys()}

# plot results
keys = ['0.95', '0.75', '0.50','0.25', '0.05']
cols = ['#CF0000', '#90FB90','#800000', '#006400']
# plot the delta EF communities with same f
fig, ax, ind = bars(EF_data, keys, col =cols[:2])
# add the deltaEF of the different f communities
fig, ax, ind = bars(EF_data_dif, keys, fig, ax, col = cols[2:])

# adjust axis
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys)
ax.set_xlabel("Proportion p of species present at both sites", fontsize = 16)
ax.set_ylim([-80,140]) # add enough space for the legend

# add legend            
legend = {col: ax.bar(0,0,width = 0 ,color = 'white'
                      ,edgecolor=col,linewidth=1.5) for col in cols}
lab = ['e>0, same f', 'e<0, same f', 'e>0, diff f', 'e<0, diff f']
ax.legend([legend[col] for col in cols],lab, loc = 'upper left')

# save figure
fig.savefig("Figure S5, diff f.pdf")

