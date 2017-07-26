"""
@author: J.W. Spaak
Contains a non vectorized version of rand_par an coex_test

rand_par randomly generates ONE community and checks whether this community
is able to coexist in the ref and the changed site (does so via coex_test)
"""

import numpy as np

from numpy.random import uniform as uni

sqrt = np.sqrt(3)

def rand_par(ave_max = 0.5, e_max = 1,
                  ave_min = -0.5,e_min = -1):
    """returns randomized parameters for one community
    
    ave_max and ave_min are the maximum/minimum that are allowed for the
    average sensitivity of each species type
    e_max, e_min are the maximum/minimum that are allowed for the sensitivity
    of each species (individually)    
    
    av means average, t means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data 6"""
    
    ave = uni(ave_min,ave_max) # average sensitivity
    alpha = uni(-0.95,-0.05) #interaction coefficient
    t_f = uni(-1/sqrt, 1/sqrt) # stdv/mean of per capita contribution
    parameters = [1,0,0,0,0,0,0] #values that do not coexist
    # runs until community that coexists is found
    while not (coex_test(*parameters)):
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        
        n = np.random.randint(5,21) # number of species
        comp = -alpha*n/(1-alpha*(n-1)) #effective competition
        t_mu = uni(comp-1,1-comp)/sqrt
        minimum = min(np.abs([e_max/ave-1, 1-e_min/ave]))
        #ensures, that sensitivity e_i  is in [e_min, e_max]
        t_e = uni(-minimum/sqrt,minimum/sqrt)
        parameters = ave,t_e,t_mu,t_f,comp,alpha,n
    return ave,t_e,t_mu,t_f,comp,alpha,n

def coex_test(ave,t_e,t_mu,t_f,comp,alpha,n):
    """tests if coexistence is given in changed site
    
    returns True, iff coexistence is guaranteed, see supplemental Info 11"""
    #computes the growthrate of species in the changed site    
    mu_change = lambda tmu, ave, te, u_i: \
                    (1+u_i*tmu*sqrt)*(1/ave-(1+u_i*te*sqrt))
    #miminum always occurs on boundaries
    min1 = mu_change(t_mu, ave, t_e,1)/(1/ave-1-t_mu*t_e)
    min2 = mu_change(t_mu, ave, t_e,-1)/(1/ave-1-t_mu*t_e)
    return min(min1, min2)>comp