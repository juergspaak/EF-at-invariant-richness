# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:46:27 2016

@author: spaakjue
Computes the DeltaEF on random parameters, trying to average them
"""

from timeit import default_timer as timer
import numpy as np

import community_construction as community

iterations = 100
save = np.zeros(iterations)
print("estimated time" , iterations*6/1000)
start = timer()
numb_par = 19
para = np.zeros([iterations,numb_par])
for i in range(iterations):
    para[i] = community.rand_par_repl()
    save[i]=min(community.rel_delta_EF_repl(*para[i]),1)
    
end=timer()
print(end-start)

paras = numb_par*[0]
for i in range(numb_par):
    paras[i] = [para[j][i] for j in range(iterations)]
labels = [r'$\bar{\mu^c}$',r'$t_{\mu^c}$', r'$\bar{\mu^r}$',r'$t_{\mu^r}$', r'$\bar{\mu^s}$',r'$t_{\mu^s}$',\
          r'$\bar{e^c}$',r'$t_{e^c}$', r'$\bar{e^s}$',r'$t_{e^s}$',r'$C_\alpha^n$',r'$t_{f^c}$',r'$t_{f^r}$',\
          r'$t_{f^s}$',"avfc","avfr","avfs",r'$\alpha$',r'$p$'] 

