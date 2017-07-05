"""
@author: J.W. Spaak
This programm plots Fig. S4
"""
import percentiles
import community_construction_repl as repl
import community_construction_coex as coex
  
keys = ['0.95', '0.75', '0.50','0.25', '0.05']

#compute DeltaEF/EF for the different communities
EF_data = {}
EF_data_dif = {}
for k in keys:
    EF_data["e<0,"+k] = repl.delta_EF_lin(*repl.para["e<0,"+k])
    EF_data["e>0,"+k] = repl.delta_EF_lin(*repl.para["e>0,"+k])
    EF_data_dif["e<0,"+k] = repl.delta_EF_lin(*repl.para["e<0,"+k], False)
    EF_data_dif["e>0,"+k] = repl.delta_EF_lin(*repl.para["e>0,"+k], False)
        

# plot results
ticks = [-100,-80,-60,-40,-20,0,20,40,60,80,100]

cols = ['#CF0000', '#90FB90','#800000', '#006400']
fig, ax, ind = percentiles.bars(EF_data, keys, col =cols[:2])
fig, ax, ind = percentiles.bars(EF_data_dif, keys, fig, ax, col = cols[2:])
ax.set_xlim([-0.2,ind[-1]+0.7])

ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys)
legend = {col: ax.bar(0,0,width = 0 ,color = 'white'
                      ,edgecolor=col,linewidth=1.5) for col in cols}
lab = ['e>0, same f', 'e<0, same f', 'e>0, diff f', 'e<0, diff f']
ax.set_xlabel("Percent of species present at both sites", fontsize = 16)
ax.set_ylim([-80,140])
ax.legend([legend[col] for col in cols],lab, loc = 'upper left')
fig.savefig("Figure newS5, diff f, 5barplots.pdf")

