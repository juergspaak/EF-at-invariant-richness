"""
@author: J.W. Spaak
This program was used to make Fig.3
"""
from percentiles import percentiles
import community_construction_repl as repl
import community_construction_coex as coex
    
keys = ['0.95', '0.50', 'rand', '0.05']
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
labels = ["p = 1.00,", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
color = ['blue','green','red', 'cyan', 'purple']
#plottig percentile curves
fig,ax = percentiles(EF_data, keys,y_min = -80, y_max = 100, 
        ticks = ticks,color = color, labels = labels)
#save figure
fig.savefig("Figure 3, changing community structure.pdf")