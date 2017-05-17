"""
@author: J.W. Spaak
This program was used to make Fig.2 compute tables S3 and S5
"""
import numpy as np
import community_construction_coex as coex
import community_construction_repl as repl

# get the distributions of the communities in invariant structure
coex_arr = {'ave': [par[0] for par in coex.para['e>0']],
                    't_e': [par[1] for par in coex.para['e>0']],
                    't_mu': [par[2] for par in coex.para['e>0']],
                    'alpha': [par[-2] for par in coex.para['e>0']],
                    'n': [par[-1] for par in coex.para['e>0']]}

temp_keys = sorted(coex_arr.keys())
coex_keys = []
for i in range(len(temp_keys)): #compute each key pair only once, cov is symmetric
    for j in range(i):
        coex_keys.append([temp_keys[i],temp_keys[j]])
tableS3 = {key1+','+key2: np.corrcoef(coex_arr[key1], coex_arr[key2])[0,1]
                    for (key1, key2) in coex_keys}
#compare to uniform distribution
com_uniform = lambda dist: np.sqrt(12)*np.std(dist)/(max(dist)-min(dist))
#insert diagional
for key in coex_arr.keys():
    tableS3[key+','+key] = com_uniform(coex_arr[key])

#get all the distributions    
rand_par = repl.para['e>0,rand']
mu = [par[0] for par in rand_par]
keys = ['avc', 'avb', 'tc' , 'tb','avu','tu']
repl_arr = {'mu,'+key: [par[key] for par in mu] for key in keys}
e = [par[1] for par in rand_par]
repl_arr.update({'e,'+key: [par[key] for par in e] for key in keys[:4]})
repl_arr['p'] = [par[-1] for par in rand_par]
repl_arr['alpha']  = [par[-2] for par in rand_par]

#compute each key pair only once
temp_keys = sorted(repl_arr.keys())
repl_keys = []
for i in range(len(temp_keys)):
    for j in range(i):
        repl_keys.append([temp_keys[i],temp_keys[j]])
#compute correlations
tableS5 = {key1+','+key2: np.corrcoef(repl_arr[key1], repl_arr[key2])[0,1]
                    for (key1, key2) in repl_keys}
#compute diagonals                   
for key in repl_arr.keys():
    tableS5[key+','+key] = com_uniform(repl_arr[key])