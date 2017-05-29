"""
@author: J.W. Spaak
This files plots Fig. S5,S6
"""
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import linregress

from percentiles import percentiles
import community_construction_coex as coex
import community_construction_repl as repl

sqrt = np.sqrt(3)
# computes the relative change in average growthrate (mu-mu')/mu
eff_ave = {'e>0,coex': coex.para['e>0'][0],
        'e<0,coex': coex.para['e<0'][0]}
keys = ['e<0,0.95', 'e<0,0.50', 'e<0,rand', 'e<0,0.05',
        'e>0,0.95', 'e>0,0.50', 'e>0,rand', 'e>0,0.05']

#same for species replacing
eff_ave_prov = np.zeros(len(repl.para[keys[0]]))
for key in keys:
    mu, e,comp,f,p, alpha = repl.para[key]
    mu_av = p*mu['avb']+(1-p)*mu['avu'] #average growth rate in ref site
    eff_ave[key] = (mu_av-mu['av_change'])/mu_av
    
keys = ['coex','0.95', '0.50', 'rand', '0.05']
labels = ["p = 1.00,", "p = 0.95", "p = 0.50", "p is random", "p = 0.05"]
color = ['blue','green','red', 'cyan', 'purple']
#plot S4
S5,ax = percentiles(eff_ave, keys, color, labels = labels)
ax.set_ylabel(r"$(\overline{\mu}-\overline{\mu'})/\overline{\mu}$",
               fontsize = 16)

key = 'rand'
adapted_para = {}
num = len(eff_ave['e<0,'+key]) #number of communities
#generate invariant community with similar effects on density by change
for change in ['e<0,', 'e>0,']:
    #linearly approximate the change in \mu
    slope, intercept, a,b,c = \
        linregress(np.linspace(0,100, num),sorted(eff_ave[change+key]))
    ave_max = 100*slope+intercept # maximum for average
    ave_min = intercept #minimum for average
    dist = min(ave_min -(-1), 1-ave_max) #distance to boundaray [-1,1]
    e_min, e_max = ave_min-dist, ave_max+dist #to have symmetric conditions
    adapted_para[change] = coex.rand_par(ave_min = ave_min, ave_max = ave_max,
                           e_min = e_min, e_max = e_max, num=num)

# in growth rates of new comunities
eff_ave.update({'e>0,ad': adapted_para['e>0,'][0],
        'e<0,ad': adapted_para['e<0,'][0]})

#compute EF for comparison
EF_data = {key+'ad': coex.delta_EF_lin(*adapted_para[key]) 
            for key in adapted_para.keys()}
EF_data.update({key: repl.delta_EF_lin(*repl.para[key]) 
            for key in repl.para.keys()})
EF_data.update({key+',coex': coex.delta_EF_lin(*coex.para[key]) 
            for key in coex.para.keys()})

#plotting
S6, (ax1,ax2) = plt.subplots(1,2,figsize =(18,7))
ref_col = color[keys.index(key)] #to have matching colors
color = [ref_col, 'orange', 'blue']
labels = ['Changing community', r'Adapted e', r'Normal e']
ticks = [-80,-60,-40,-20,0,20,40,60,80,100]
percentiles(eff_ave, [key, 'ad', 'coex'], color,labels = labels,
                y_max = 0.6, ticks = ticks, plot = (S6, ax1))
ticks = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
percentiles(EF_data, [key, 'ad', 'coex'], color,labels = labels,
                y_max = 100, plot = (S6, ax2))
ax1.set_xlabel("percentile", fontsize = 16)
ax1.set_ylabel(r"$(\overline{\mu}-\overline{\mu'})/\overline{\mu}$",
               fontsize = 16)
ax1.set_title('A. Change in growthrates', fontsize = 16)
ax2.set_title('B. Ecosystem functioning', fontsize = 16)
ax2.set_ylabel(r'$100\cdot\Delta EF/EF$', fontsize = 16)
#save figure
S5.savefig("Figure S5, Adapted e.pdf")
S6.savefig("Figure S6, Adapted e.pdf")

    
