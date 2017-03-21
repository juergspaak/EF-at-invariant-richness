# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:07:05 2017

@author: spaakjue
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 12:25:08 2017

@author: spaakjue
"""
npz = np.load("npz, save_all repl.npz")
coex = np.load("npz, species_coex.npz")
npz_avmue = np.load("npz, adapted mue-mu.npz")
avmue_div_mu={}
keys = ['paras_p05',  'paras_p50',  'paras_p95',
       'paras_prand']
deltl = {}
for key in keys:
    avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p  = npz[key]
    a = p*avmuc*avec*(1/avec-1-t_muc*t_ec)+(1-p)*avmus*aves*(1/aves-1-t_es*t_mus)
    b = p*avmuc +(1-p)*avmur
    avmue_div_mu[key] = a/b
    deltl[key] = 100*npz['delta_EF_p'+key[7:]]

ave,t_e,t_mu,t_f,comp,alpha,n = coex['paras']
a = ave*(1/ave-1-t_e*t_mu)
avmue_div_mu['coex'] = a
deltl['coex'] = 100*npz['delta_EF_p100']

key_ord = ['coex', 'paras_p95', 'paras_p50', 'paras_prand', 'paras_p05']
yl  = r'$\overline{\mu(1-e)}/\overline{\mu}$'
lab = ["p = 1.00", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
plot_percentiles(avmue_div_mu,  y_label = yl,ordered_keys = key_ord, grid = False,
                 labels = lab, save = "avmue_dov_avmu2")
yl = r'$100\cdot\Delta EF/EF$'
lab2 = [r'adapted $max(\bar{e})\approx 0.55$', 'p = 0.50']

plot_percentiles([npz_avmue['delta_EF'], 100*npz['delta_EF_p50']],  y_label = yl, grid = False,
                 labels = lab2, save = "adapted avemax2",color = ['b','r'])

plot_percentiles(deltl, labels = lab, y_label = yl, grid = False, save = "delta_EF, different p,2", ordered_keys = key_ord)


""" for the covariances"""
def fun_covs(para, delta_EF):
    covs = {'11':[], '1-1':[], '-11':[], '-1-1':[]}
    for i in range(len(delta_EF)):
        a = str(int(np.sign(para[i][1]*para[i][2])))
        b = str(int(np.sign(para[i][1]*para[i][3])))
        covs[a+b].append(delta_EF[i])
    covs2 = {}
    covs2[r'cov$(e,\mu) = 1$, cov$(e,f)=1$'] = covs['11']
    covs2[r'cov$(e,\mu) = -1$, cov$(e,f)=-1$'] = covs['-1-1']
    covs2[r'cov$(e,\mu) = -1$, cov$(e,f)=1$'] = covs['-11']
    covs2[r'cov$(e,\mu) = 1$, cov$(e,f)=-1$'] = covs['1-1']
    return covs2
    
covs = fun_covs(coex['para'], coex['delta_EF_coex'])
k2,k3, k1, k4= sorted(covs.keys())
ord_key = [k1,k2,k3,k4]
ls = ['-', '--', '-.',':']
plot_percentiles(covs, ordered_keys = ord_key, y_label = yl, grid = False, save = "covariances2",
                 ls = ls)

n = 20
""" for the adjustment terms:"""
adjust = {}
for p in ['05','50','95','rand']:
    adjust[p] = []
    for i in npz['para_p'+p]:
        adjust[p].append(rel_delta_EF_repl(*i,0))

color = {'95':['g','g'],'50':['r','r'],'rand':['c','c'],'05':['m','m']} 
ls = [':','--']
for p in ['05','50','95','rand']:
    labels = ["without adjustment terms, p = 0."+p, "with adjustment terms, p = 0."+p]
    if p == 'rand':
        labels = ["without adjustment terms, p is random", "with adjustment terms, p is random"]
    plot_percentiles([adjust[p], 100*npz['delta_EF_p'+p]], y_label = yl,
                     y_max = 100, save = "adjust,"+p, grid = False,
                     labels = labels, ls = ls,color = color[p])

""" for the distribution of the parameters in all coex"""
plot_percentiles([coex['paras'][3], coex['paras'][1]], y_label = r'$t_e$', grid = False,
                 labels = ["uniform distribution", "effective distribution"],save = "t_e")
plot_percentiles([coex['paras'][3], coex['paras'][2]], y_label = r'$t_{\mu}$', grid = False,
                 labels = ["uniform distribution", "effective distribution"], save = "t_mu")
plot_percentiles([np.random.randint(5,21,int(1e6)), coex['paras'][-1]], y_label = r'$n$',
                 grid = False, text = "hi",
                 labels = ["uniform distribution", "effective distribution"],save ="n")
    
