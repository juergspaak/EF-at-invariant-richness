"""
@author: J.W. Spaak
This programm plots Fig. 4
"""
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

# plot results
keys = ['ref'] + keys
fig, ax, ind = percentiles.bars(EF_data, keys)

# adjust axis
ax.set_xlim([-0.2,ind[-1]+0.7])
ax.set_xlabel(r"Maximal value for $\bar{H}$", fontsize = 16)
ax.set_xticks(ind+0.15)
ax.set_xticklabels(keys)
ax.set_ylim(-60, 80)

# save the figure
fig.savefig("Figure S4, asymptotic fucntioning.pdf")