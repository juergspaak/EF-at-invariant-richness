"""
@author: J.W. Spaak
This programm plots Fig. 2
"""
import numpy as np

from percentiles import percentiles
import community_construction_coex as coex
    
#split into different covariances
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
                
EF_covs = {key: coex.delta_EF_lin(*(par_covs[key])) for key in par_covs.keys()}       
# plot results
ticks = [-60,-40,-20,0,20,40,60,80,100]

labels = [r'cov$(|e|,\mu) = 1$, cov$(e,f)=1$',
              r'cov$(|e|,\mu) = -1$, cov$(e,f)=1$',
            r'cov$(|e|,\mu) = -1$, cov$(e,f)=-1$',
            r'cov$(|e|,\mu) = 1$, cov$(e,f)=-1$']
color = ['purple','red', 'green', 'blue']
#plot percentile curves
fig,ax = percentiles(EF_covs, keys, y_max = 100, ticks = ticks,
        color = color, labels = labels)
#save figure
fig.savefig("Figure 2, Covariances.pdf")
