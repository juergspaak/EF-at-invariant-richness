"""
@author: J.W. Spaak
Contains all functions, that are needed to generate communities as well as
compute the EF of those communities

rand_par randomly generates a community and checks whether this community
is able to coexist in the ref and the changed site(does so via coex_test)

delta_EF_lin computes \DeltaEF/EF with linear contribution to function
delta_EF_asym with asymptotic function
"""

import numpy as np

from numpy.random import uniform as uni
from scipy.integrate import quad

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
    parameters = [1,0,0,0,0,0,0] #values that do not coexist
    # runs until community that coexists is found
    # WARNING: changing settings might turn this into an infinite loop
    while not (coex_test(*parameters)):
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        t_f = uni(-1/sqrt, 1/sqrt) # stdv/mean of per capita contribution
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
    
def delta_EF_lin(ave,t_e,t_mu,t_f,comp,alpha,n):
    """ computes DetalEF/EF for one community
    
    for computational background see Eq. 4"""
    save_1 = -(1+t_mu*t_e)*ave
    save_2 = t_f*(t_mu+t_e)/(1+t_mu*t_e)-t_f*t_mu
    save_3 = t_mu*t_f+1-comp
    return 100*save_1*(1+save_2/save_3)
    
def delta_EF_asym(ave,t_e,t_mu,t_f,comp,alpha,n,max_ave_H=1):
    """computes the EF with asymptotic f, f(N) = f_i*H_i*N_i/(N_i+H_i)
    
    ave,t_e,t_mu,t_f,comp,alpha,n should be the return values of rand_par_coex
    max_ave_H is the maximum value for average value for  H.
    H_i is uniformly distributed in [0,2*ave_H]"""
    # choose distribution of H: H ~u[0,2*ave]
    ave_H = uni(0,max_ave_H)
    t_H = uni(-1/sqrt, 1/sqrt)
    H = lambda x: ave_H*(1+t_H*sqrt*x)
    #asymptotic EF in N, f(N) = f_i*H_i*N_i/(N_i+H_i)
    eco_fun = lambda x, N: n*(1+t_f*x*sqrt)*H(x)*N(x)/(N(x)+H(x))
    # computes the equilibrium densities of species N, in changed and ref site
    N_ref = lambda x: (1+t_mu*sqrt*x-comp)/(1+alpha)
    N_change = lambda x: ((1+x*t_mu*sqrt)*(1-ave*(1+t_e*sqrt*x))-\
                comp*(1-ave*(1+t_mu*t_e)))/(1+alpha)
    
    # Integrate for the average, divide by 2 for normalisation
    EF_ref = quad(lambda x: eco_fun(x,N_ref)/2,-1,1)[0] 
    EF_change = quad(lambda x: eco_fun(x,N_change)/2,-1,1)[0]
    return 100*(EF_change-EF_ref)/EF_ref #multiply by 100 for percent
