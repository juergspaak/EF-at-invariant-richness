"""
@author: J.W. Spaak
This program computes tables S3 and S5
"""
import numpy as np
from pandas import DataFrame
import community_construction_coex as coex
import community_construction_repl as repl

# get the distributions of the communities in invariant structure
coex_arr = {'ave': coex.para['e>0'][0],
            't_e': coex.para['e>0'][1],
            't_mu': coex.para['e>0'][2],
            'alpha': coex.para['e>0'][-1],
            'n': coex.para['e>0'][-2]}

temp_keys = sorted(coex_arr.keys())
coex_keys = []
for i in range(len(temp_keys)): 
    for j in range(len(temp_keys)):
        coex_keys.append([temp_keys[i],temp_keys[j]])
valuesS3 = {key1+','+key2: np.corrcoef(coex_arr[key1], coex_arr[key2])[0,1]
                    for (key1, key2) in coex_keys}
#compare to uniform distribution
com_uniform = lambda dist: np.sqrt(12)*np.std(dist)/(max(dist)-min(dist))
#insert diagional
for key in coex_arr.keys():
    valuesS3[key+','+key] = com_uniform(coex_arr[key])
#printing the table S3
ord_keys = ["ave", "t_e", "t_mu", "alpha", "n"]
matS3 = np.zeros([len(ord_keys), len(ord_keys)])
for i in range(len(ord_keys)):
    for j in range(len(ord_keys)):
        if np.abs(valuesS3[ord_keys[i]+','+ord_keys[j]])>0.05:
            matS3[i,j] = valuesS3[ord_keys[i]+','+ord_keys[j]]
        
tableS3 = DataFrame(np.round(matS3,2), ord_keys, ord_keys)
print("Table S3")
print(tableS3)
#get all the distributions    
mu, e,comp,f, p, alpha  = repl.para['e>0,0.50']
keys = ['avc', 'avb', 'tc' , 'tb','avu','tu']
repl_arr = {'mu,'+key: mu[key] for key in keys}
repl_arr.update({'e,'+key: e[key] for key in keys[:4]})
repl_arr.update({"alpha": alpha, "p":p})

#compute each key pair only once
temp_keys = sorted(repl_arr.keys())
repl_keys = []
for i in range(len(temp_keys)):
    for j in range(len(temp_keys)):
        repl_keys.append([temp_keys[i],temp_keys[j]])
#compute correlations
valuesS5 = {key1+','+key2: np.corrcoef(repl_arr[key1], repl_arr[key2])[0,1]
                    for (key1, key2) in repl_keys}
#compute diagonals                   
for key in repl_arr.keys():
    valuesS5[key+','+key] = com_uniform(repl_arr[key])
    
for key in coex_arr.keys():
    valuesS3[key+','+key] = com_uniform(coex_arr[key])
#printing the table S5    
ord_keys = temp_keys[5:-2]+temp_keys[1:5]+['alpha','p']
matS5 = np.zeros([len(ord_keys), len(ord_keys)])
for i in range(len(ord_keys)):
    for j in range(len(ord_keys)):
        if np.abs(valuesS5[ord_keys[i]+','+ord_keys[j]])>0.05:
            matS5[i,j] = valuesS5[ord_keys[i]+','+ord_keys[j]]
    
tableS5 = DataFrame(np.round(matS5,2), ord_keys, ord_keys)
print("\nTable S5")
print(tableS5.dropna())
    
