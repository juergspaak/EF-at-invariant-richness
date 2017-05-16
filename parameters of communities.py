"""
@author: J.W. Spaak
This files generates the communities and saves them
"""
import pickle

import help_functions as hef
import community_construction_repl as repl
import community_construction_coex as coex

#generate communities, invariant species composition
num_com = int(1e5) #number of communities computed
com_para_coex = {}
com_para_coex['e>0'] = hef.com_con(coex.rand_par,num_com,ave_min=0,e_min=0)
com_para_coex['e<0'] = hef.com_con(coex.rand_par,num_com,ave_max=0,e_max=0)
#save data for future use
pickle.dump(com_para_coex, open("coex, com_para.p","wb"))

#generate communities for different p
keys = ['0.95', '0.50', 'rand', '0.05']
com_para_repl = {}
for key in keys:
    try:
        p= float(key)
    except ValueError:
        p = key
    com_para_repl['e>0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_min=0,e_min=0)
    com_para_repl['e<0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_max=0,e_max=0)
#save data for future use  
pickle.dump(com_para_repl, open("repl, com_para.p","wb"))
