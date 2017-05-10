# -*- coding: utf-8 -*-
"""
Created on Thu Dec 15 09:46:27 2016

@author: spaakjue
Computes the DeltaEF on random parameters, trying to average them
"""

from timeit import default_timer as timer
import numpy as np

import community_construction as community
def run_programm(p, iterations= int(1e3), savenpz =  False):
    save = np.zeros(iterations)
    print("estimated time" , iterations*6/1000)
    start = timer()
    numb_par = 19
    para = np.zeros([iterations,numb_par])
    for i in range(iterations):
        para[i] = community.rand_par_repl(p = p, ave_max = 0, e_max = 0)
        save[i]=community.rel_delta_EF_repl(*para[i])
        if (i % 10000)==0:
            print(p,i, "done")
            if savenpz:
                np.savez("npz, repl, neg e, p="+str(p)+",savety.npz",
                     rel_delta_EF = save[:i], para = para[:i])
        
    end=timer()
    print("time needed for p="+str(p)+","+str(end-start), "\n")
    
    paras = numb_par*[0]
    for i in range(numb_par):
        paras[i] = [para[j][i] for j in range(iterations)]
    labels = [r'$\bar{\mu^c}$',r'$t_{\mu^c}$', r'$\bar{\mu^r}$',r'$t_{\mu^r}$', r'$\bar{\mu^s}$',r'$t_{\mu^s}$',\
              r'$\bar{e^c}$',r'$t_{e^c}$', r'$\bar{e^s}$',r'$t_{e^s}$',r'$C_\alpha^n$',r'$t_{f^c}$',r'$t_{f^r}$',\
              r'$t_{f^s}$',"avfc","avfr","avfs",r'$\alpha$',r'$p$']
    if savenpz:
        np.savez("npz, repl, neg e, p="+str(p)+".npz",
                    paras =paras,rel_delta_EF = save, para = para)
    return save,para, paras
    
    
save05,para05,paras05 = run_programm(0.05)
save95,para95,paras95 = run_programm(0.95)
save50 ,para50,paras50= run_programm(0.5)
saverand ,pararand,parasrand= run_programm('rand')
labels = [r'$\bar{\mu^c}$',r'$t_{\mu^c}$', r'$\bar{\mu^r}$',r'$t_{\mu^r}$', r'$\bar{\mu^s}$',r'$t_{\mu^s}$',\
              r'$\bar{e^c}$',r'$t_{e^c}$', r'$\bar{e^s}$',r'$t_{e^s}$',r'$C_\alpha^n$',r'$t_{f^c}$',r'$t_{f^r}$',\
              r'$t_{f^s}$',"avfc","avfr","avfs",r'$\alpha$',r'$p$']
plot_percentiles([save05,save50,save95,saverand], [5,50,95,"rand"])
"""
npz05 = np.load("npz, repl, neg e, p=05.npz")
npz50 = np.load("npz, repl, neg e, p=50.npz")
npz95 = np.load("npz, repl, neg e, p=95.npz")
npzrand = np.load("npz, repl, neg e, p=rand.npz")
npz100 = np.load("npz, repl, neg e, p=100.npz")
np.savez("npz, save_all repl, neg ave.npz",
         delta_EF_p05 = npz05['rel_delta_EF'], paras_p05 = npz05['paras'],
        para_p05 = npz05['para'],
        delta_EF_p50 = npz50['rel_delta_EF'], paras_p50 = npz50['paras'],
        para_p50 = npz50['para'],
        delta_EF_p95 = npz95['rel_delta_EF'], paras_p95 = npz95['paras'],
        para_p95 = npz95['para'],
        delta_EF_prand = npzrand['rel_delta_EF'], paras_prand = npzrand['paras'],
        para_prand = npzrand['para'],
        delta_EF_p100 = npz100['rel_delta_EF'], paras_p100 = npz100['paras'],
        para_p100 = npz100['para'])"""

