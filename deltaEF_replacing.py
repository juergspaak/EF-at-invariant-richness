# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:46:27 2016

@author: spaakjue
Computes the DeltaEF on random parameters, trying to average them
"""

from timeit import default_timer as timer
import numpy as np

import community_construction as community
def run_programm(p, iterations= int(1e3)):
    save = np.zeros(iterations)
    print("estimated time" , iterations*6/1000)
    start = timer()
    numb_par = 19
    para = np.zeros([iterations,numb_par])
    for i in range(iterations):
        para[i] = community.rand_par_repl(p = p, ave_max = 0, e_max = 0)
        save[i]=community.rel_delta_EF_repl(*para[i])
        if (i % 10000)==0:
            print("\n", p,i, "done\n")
            #np.savez("npz2, repl, neg e, p="+str(p)+",savety.npz",
            #         rel_delta_EF = save[:i], para = para[:i])
        
    end=timer()
    print("time needed for p="+str(p)+","+str(end-start))
    
    paras = numb_par*[0]
    for i in range(numb_par):
        paras[i] = [para[j][i] for j in range(iterations)]
    labels = [r'$\bar{\mu^c}$',r'$t_{\mu^c}$', r'$\bar{\mu^r}$',r'$t_{\mu^r}$', r'$\bar{\mu^s}$',r'$t_{\mu^s}$',\
              r'$\bar{e^c}$',r'$t_{e^c}$', r'$\bar{e^s}$',r'$t_{e^s}$',r'$C_\alpha^n$',r'$t_{f^c}$',r'$t_{f^r}$',\
              r'$t_{f^s}$',"avfc","avfr","avfs",r'$\alpha$',r'$p$']
    #np.savez("npz2, repl, neg e, p="+str(p)+".npz",
    #                 paras =paras,rel_delta_EF = save, para = para)
    return save,para, paras
    
    
save05,para,paras = run_programm(0.05)
save95,para,paras = run_programm(0.95)
save50 ,para,paras= run_programm(0.5)
saverand ,para,paras= run_programm('rand')
labels = [r'$\bar{\mu^c}$',r'$t_{\mu^c}$', r'$\bar{\mu^r}$',r'$t_{\mu^r}$', r'$\bar{\mu^s}$',r'$t_{\mu^s}$',\
              r'$\bar{e^c}$',r'$t_{e^c}$', r'$\bar{e^s}$',r'$t_{e^s}$',r'$C_\alpha^n$',r'$t_{f^c}$',r'$t_{f^r}$',\
              r'$t_{f^s}$',"avfc","avfr","avfs",r'$\alpha$',r'$p$']
