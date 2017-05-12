"""
Created on Wed Feb  8 13:07:05 2017

@author: spaakjue
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 12:25:08 2017

@author: spaakjue
"""
try:
    loaded==True
except NameError:
    npz = np.load("npz, save_all repl.npz")
    npz_neg = np.load("npz, save_all repl, neg ave.npz")
    coex = np.load("npz, species_coex.npz")
    npz_avmue = np.load("npz, adapted mue-mu.npz")
    loaded = True
    avmue_div_mu={}
    keys = ['paras_p05',  'paras_p50',  'paras_p95',
       'paras_prand']
    deltl = {'05neg': npz_neg['delta_EF_p05'], '50neg': npz_neg['delta_EF_p50'],
          '95neg': npz_neg['delta_EF_p95'],'100neg': npz_neg['delta_EF_p100']
            ,'randneg': npz_neg['delta_EF_prand']}
    for key in keys:
        avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
                  t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
                  avfc,avfr,avfs,alpha,p  = npz[key]
        a = p*avmuc*avec*(1/avec-1-t_muc*t_ec)+(1-p)*avmus*aves*(1/aves-1-t_es*t_mus)
        b = p*avmuc +(1-p)*avmur
        avmue_div_mu[key] = a/b
        deltl[key] = 100*npz['delta_EF_p'+key[7:]]
    yl = r'$100\cdot\Delta EF/EF$'
    ave,t_e,t_mu,t_f,comp,alpha,n = coex['paras']
    a = ave*(1/ave-1-t_e*t_mu)
    avmue_div_mu['coex'] = a
    deltl['coex'] = 100*npz['delta_EF_p100']
    """ for the covariances"""
    def fun_covs(para, delta_EF, adinfo):
        covs = {'11':[], '1-1':[], '-11':[], '-1-1':[]}
        for i in range(len(delta_EF)):
            a = str(int(np.sign(para[i][1]*para[i][2])))
            b = str(int(np.sign(para[i][1]*para[i][3])))
            covs[a+b].append(delta_EF[i])
        covs2 = {adinfo+key: covs[key] for key in covs.keys()}
        return covs2
    covs = fun_covs(coex['para'], coex['delta_EF_coex'], "e>0")
    covs.update(fun_covs(npz_neg['para_p100'], npz_neg['delta_EF_p100'],"e<0"))
    marker = ['^','^','^','^','^','v','v','v','v','v']

def avmue_dov_avmu2():
    key_ord = ['coex', 'paras_p95', 'paras_p50', 'paras_prand', 'paras_p05']
    yl  = r'$\overline{\mu(1-e)}/\overline{\mu}$'
    lab = ["p = 1.00,", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
    plot_percentiles(avmue_div_mu,  y_label = yl,ordered_keys = key_ord, grid = False,
                     labels = lab, save = "avmue_dov_avmu2")
def adapted_avemax2():
    lab2 = [r'adapted $max(\bar{e})\approx 0.55$', 'p = 0.50']
    key_ord = ['coex', 'paras_p95', 'paras_p50', 'paras_prand', 'paras_p05',
               '100neg', '95neg', '50neg', 'randneg', '05neg']
    
    plot_percentiles([npz_avmue['delta_EF'], 100*npz['delta_EF_p50']],  y_label = yl, grid = False,
                     labels = lab2, save = "adapted avemax2",color = ['b','r'])

def deltl_plot():
    color = 2*["blue", "green", "red", "cyan", "magenta"]
    k = sorted(deltl.keys())
    ord_keys = [k[1], k[3], k[2], k[9],k[0], k[4],k[7], k[6], k[8], k[5]]
    fig,ax = plot_percentiles(deltl, y_label = yl, grid = False, 
                     y_max = 100,marker =marker,color = color,
                     ordered_keys = ord_keys)
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='^', linestyle='')
    anyArtist = plt.Line2D((0,2),(0,0), color='k',marker='v', linestyle='')
    handles, labels = ax.get_legend_handles_labels()
    display = range(5)
    for i in range(5):
        handles[i] = plt.Line2D((0,2),(0,0), color=color[i])
    labels = ["p = 1.00,", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
    ax.legend(loc = "upper left")
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[simArtist,anyArtist],
              labels+['e<0', 'e>0']
               ,loc = "upper left", numpoints = 1)
    fig.savefig("Figure, percentiles delta_EF, different p,2.pdf", dpi=fig.dpi)
    
def covs_plot(structure = "reversed"):
    colors = 2*['purple', 'red', 'green', 'darkblue']
    k= sorted(covs.keys())
    if structure == "pointsymmetry":
        ord_key = [k[1], k[3],k[2],k[0],k[7],k[5],k[4],k[6]]
        labels = [r'cov$(|e|,\mu) = 1$, cov$(e,f)=1$',
              r'cov$(|e|,\mu) = -1$, cov$(e,f)=1$',
            r'cov$(|e|,\mu) = -1$, cov$(e,f)=-1$',
            r'cov$(|e|,\mu) = 1$, cov$(e,f)=-1$']
    elif structure == "reversed":
        ord_key = [k[1], k[3],k[2],k[0],k[6],k[4],k[5],k[7]]
        labels = [r'cov$(|e|,\mu) = 1$, cov$(|e|,f)=-1$',
              r'cov$(|e|,\mu) = -1$, cov$(|e|,f)=-1$',
            r'cov$(|e|,\mu) = -1$, cov$(|e|,f)=1$',
            r'cov$(|e|,\mu) = 1$, cov$(|e|,f)=1$']
    else:
        ord_key = [k[1], k[3],k[2],k[0],k[5],k[7],k[6],k[4]]
        labels = [r'cov$(e,\mu) = -1$, cov$(e,f)=1$',
              r'cov$(e,\mu) = 1$, cov$(e,f)=1$',
            r'cov$(e,\mu) = 1$, cov$(e,f)=-1$',
            r'cov$(e,\mu) = -1$, cov$(e,f)=-1$']
        structure = "normal"
    
    marker = ['^','^','^','^','v','v','v','v']
    colors = 2*['purple', 'red', 'green', 'darkblue']
    fig,ax = plot_percentiles(covs, ordered_keys = ord_key, y_label = yl,
                y_max = 100, grid = False, 
                save = "covariances2",marker =marker, color = colors)
    simArtist = plt.Line2D((0,1),(0,0), color='k', marker='^', linestyle='')
    anyArtist = plt.Line2D((0,2),(0,0), color='k',marker='v', linestyle='')
    handles, tlabels = ax.get_legend_handles_labels()
    display = range(4)
    for i in range(4):
        handles[i] = plt.Line2D((0,2),(0,0), color=colors[i])
    
    ax.legend(loc = "upper left")
    ax.legend([handle for i,handle in enumerate(handles) if i in display]+[simArtist,anyArtist],
              labels+['e<0', 'e>0']
        ,loc = "upper left", numpoints = 1)
    fig.savefig("Figure, percentiles covs, "+structure +".pdf", dpi=fig.dpi)
    
covs_plot("reversed")

n = 20
""" for the adjustment terms:"""
if (False):
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
                     grid = False,
                     labels = ["uniform distribution", "effective distribution"],save ="n")
    
