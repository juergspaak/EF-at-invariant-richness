"""
@author: J.W. Spaak
This program was used to make Fig.2
"""
import numpy as np

import help_functions as hef
import community_construction_repl as repl
import community_construction_coex as coex

num_com = int(1e4) #number of communities computed

#generate communities for different p
keys = ['0.95', '0.50', 'rand', '0.05']
com_para = {}
for key in keys:
    print(key)
    try:
        p= float(key)
    except ValueError:
        p = key
    com_para['e>0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_min=0,e_min=0)
    com_para['e<0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_max=0,e_max=0)

print("com_con")

#compute DeltaEF/EF for the different communities
EF_data = {}
for k in keys:
    EF_data["e<0,"+k] = \
        hef.comp_EF(repl.delta_EF_asym,com_para["e<0,"+k])
    EF_data["e>0,"+k] = \
        hef.comp_EF(repl.delta_EF_asym,com_para["e>0,"+k])
        
# generat referenc, with invariant community structure
com_para['e>0,ref'] = hef.com_con(coex.rand_par,num_com,ave_min=0,e_min=0)
com_para['e<0,ref'] = hef.com_con(coex.rand_par,num_com,ave_max=0,e_max=0)
# compare to linear functioning        
EF_data["e<0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,com_para["e<0,ref"])
EF_data["e>0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,com_para["e>0,ref"])
        
keys = ['ref']+keys
# plot results
ticks = [-60,-40,-20,0,20,40,60,80,100]
labels = ["p = 1.00,", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
color = ['purple','red', 'green', 'cyan','blue']
fig,ax = hef.percentiles(EF_data, keys, y_max = 100, ticks = ticks,
        color = color, labels = labels)
