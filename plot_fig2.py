"""
@author: J.W. Spaak
This program was used to make Fig.2
"""
import numpy as np

import help_functions as hef
import community_construction_coex as coex


#ensure that loading worked
try:
    coex.para
except AttributeError:
    coex.error()
    
#split into different covariances
EF_covs = {'e<0,11':[], 'e<0,1-1':[], 'e<0,-11':[], 'e<0,-1-1':[],
           'e>0,11':[], 'e>0,1-1':[], 'e>0,-11':[], 'e>0,-1-1':[]}
abs_sign = {'e>0,':1, 'e<0,':-1} #sign change for |e|
for key in ['e>0,', 'e<0,']:
    for para in coex.para[key[:3]]:
            a = str(int(np.sign(abs_sign[key]*para[1]*para[2])))
            b = str(int(np.sign(para[1]*para[3])))
            EF_covs[key+a+b].append(coex.delta_EF_lin(*para))

        
# plot results
ticks = [-60,-40,-20,0,20,40,60,80,100]
keys = ['11','-11','-1-1','1-1']
labels = [r'cov$(|e|,\mu) = 1$, cov$(e,f)=1$',
              r'cov$(|e|,\mu) = -1$, cov$(e,f)=1$',
            r'cov$(|e|,\mu) = -1$, cov$(e,f)=-1$',
            r'cov$(|e|,\mu) = 1$, cov$(e,f)=-1$']
color = ['purple','red', 'green', 'blue']
#plot percentile curves
fig,ax = hef.percentiles(EF_covs, keys, y_max = 100, ticks = ticks,
        color = color, labels = labels)
#save figure
fig.savefig("Figure 2, Covariances.pdf")
