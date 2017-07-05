"""
@author: J.W. Spaak
This programm plots Fig.S4
"""
import numpy as np
import matplotlib.pyplot as plt

#compute competition
comp=lambda n,alpha: -alpha*n/(1-alpha*(n-1))

#sample alphas
alphas = np.linspace(-1,0,1000)
fig, ax = plt.subplots()
for n in [2,5,10,15,20,100]:
    comps = [comp(n,alpha) for alpha in alphas]
    plt.plot(alphas,comps,label='n = ' + str(n))

# add legend and labels
plt.legend(loc='lower left')
ax.set_xlabel(r'$\alpha$',fontsize=14)
ax.set_ylabel(r'$C^n_{\alpha}$',fontsize=14)

#save figure
fig.savefig("Figure S1, competition.pdf")
