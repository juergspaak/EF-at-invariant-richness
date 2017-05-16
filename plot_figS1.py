# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 08:06:03 2016

@author: spaakjue
"""
import numpy as np
import matplotlib.pyplot as plt

#compute competition
comp=lambda n,alpha: -alpha*n/(1-alpha*(n-1))

#sample alphas
alphas=np.linspace(-1,0,1000)
fig=plt.figure()
for n in [2,5,10,15,20,100]:
    comps=[comp(n,alpha) for alpha in alphas]
    plt.plot(alphas,comps,label='n = ' + str(n))
  
plt.legend(loc='lower left')
fig.gca().set_xlabel(r'$\alpha$',fontsize=14)
fig.gca().set_ylabel(r'$C^n_{\alpha}$',fontsize=14)

plt.show(fig)

#save figure
fig.savefig("Figure S1, competition.pdf")
