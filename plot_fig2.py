"""
@author: J.W. Spaak
This programm plots Fig. 2
"""
import numpy as np

import percentiles
import community_construction_coex as coex
   
# split into different covariances
par_covs = {}
abs_sign = {'e>0,':1, 'e<0,':-1} #sign change for |e|
for key in ['e>0,', 'e<0,']:
    par = coex.para[key[:3]]
    a = np.sign(abs_sign[key]*par[1]*par[2]) #contains values -1,1
    b = np.sign(par[1]*par[4]) #contains values -1,1
    cases = 10*a+b # all possibilities are -11,-9,9,11
    keys = ['11','-9','-11','9']
    par_covs.update({key+case: [coex.para[key[:3]][i][cases==float(case)] 
                for i in range(7)] for case in keys})
# compute deltaEF/EF                
EF_covs = {key: coex.delta_EF_lin(*(par_covs[key])) for key in par_covs.keys()}       

# plot results
fig, ax, ind = percentiles.bars(EF_covs, keys)

# limits
ax.set_ylim(-60,100)
ax.set_xlim([-0.2,ind[-1]+0.7])

# Label at lower x-axis 
ax.set_xticks(ind+0.15)
ax.set_xticklabels(2*[r'cov$(|e|,\mu) = 1$',r'cov$(|e|,\mu) = -1$'],
                   fontsize = 14)
ax.set_xlabel("Covariances of Sensitivity and Growthrates", fontsize = 16)

# label at upper x-axis
ax2 = ax.twiny()
ax2.set_xticks([0.25,0.75])
ax2.set_xticklabels([r'cov$(|e|,f) = 1$',r'cov$(|e|,f) = -1$']
                    , fontsize = 14)
ax2.set_xlabel("Covariances of Sensitivity and per-capita Contribution", 
               fontsize = 16)

# save the figure
fig.savefig("Figure 2, Covariances, barplots.pdf")