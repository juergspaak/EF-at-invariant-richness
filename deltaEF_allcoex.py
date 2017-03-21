# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:46:27 2016

@author: spaakjue
Computes the DeltaEF on random parameters, trying to average them
The distributions try to be similar to the ones in the replacing programm
"""
import numpy as np
from timeit import default_timer as timer

import community_construction as community

counter1,counter2=0,0 #count # ecosystems with positive/negative change in EF
iterations = int(100000)   #number of random ecosystems
rel_delta_EF = np.zeros([iterations]) #saves \DeltaEF/EF
neg_para = np.zeros([iterations,7]) # saves the parameters in negative change ecosystems
neg_delta_EF = np.zeros(iterations) # saves rel_delta_EF in cases it is negative
pos_para = np.zeros([iterations,7]) # saves the parametesr in positive change ecosystems
pos_delta_EF = np.zeros(iterations) # saves rel_delta_EF in cases it's positive
para = np.zeros([iterations,7])     # saves the parameters of all ecosystems
label_names = [r'$\bar{e}$', r'$t_e$',r'$t_{\mu}$',
                r'$t_f$', r'$C_{\alpha}^n$', r'$\alpha$', r'$n$'] # names for axes used in plit function
label_delEF = r'$\Delta EF/EF$'
sqrt = np.sqrt(3)


axis = [[0,1],[0,1/np.sqrt(3)], [-1/np.sqrt(3),1/np.sqrt(3)], [-1/np.sqrt(3),1/np.sqrt(3)], 
        [0,1], [-1,0], [4.9,20.1]] #axis for plotting heatplots 

print("estimated time: ", iterations /10000*2)
start=timer()
for i in range(iterations):
    parameters = community.rand_par_coex() # generate random community
    para[i]=parameters #save all parameters
    rel_delta_EF[i] = community.rel_delta_EF_coex(*parameters)
    if rel_delta_EF[i] < -100: #mistake happend
        print(rel_delta_EF[i], parameters)
    elif rel_delta_EF[i] < 0:
        neg_delta_EF[counter1] = rel_delta_EF[i]
        neg_para[counter1,:] = para[i,:]
        counter1+=1
    else:
        pos_delta_EF[counter2] = rel_delta_EF[i]
        pos_para[counter2,:] = para[i,:]
        counter2+=1
    
        
end=timer()
print(end-start)
neg_delta_EF = neg_delta_EF[:counter1] # chop of unimportant
neg_para = neg_para[:counter1]
pos_delta_EF = pos_delta_EF[:counter2]
pos_para = pos_para[:counter2]
paras = 7*[0]   # are used to save the parameters in one list
neg_paras = 7*[0]
pos_paras = 7*[0]
for i in range(7):
    paras[i] = [para[j][i] for j in range(iterations)]
    neg_paras[i] = [neg_para[j][i] for j in range(counter1)]
    pos_paras[i] = [pos_para[j][i] for j in range(counter2)]
                                
rel_min_comp = [1-paras[4][i] for i in range(len(paras[4]))]
pos_min_comp = [1-pos_paras[4][i] for i in range(len(pos_paras[4]))]
neg_min_comp = [1-neg_paras[4][i] for i in range(len(neg_paras[4]))]

zero_delta_EF = [min(0,i) for i in rel_delta_EF]