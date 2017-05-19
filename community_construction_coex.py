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
import pickle

from numpy.random import uniform as uni
from scipy.integrate import quad

sqrt = np.sqrt(3)
#load the communities
try:
    para = pickle.load(open("coex, com_para.p", "rb"))
except FileNotFoundError: #other functions should be still be available
    def error():
        raise FileNotFoundError("No such file or directory: 'coex, com_para.p'"
            "\nPlease run the file 'parameters of communities.py "+
            "to create this file. Read the readme for further instructions.")
        
def vec_uni(low, high):
    try:
        return low+(high-low)*uni(size = len(low))
    except TypeError:
        return uni(low, high)

def rand_par(num = 1000, ave_max = 0.5, e_max = 1,
                  ave_min = -0.5,e_min = -1):
    """returns randomized parameters for one community
    
    ave_max and ave_min are the maximum/minimum that are allowed for the
    average sensitivity of each species type
    e_max, e_min are the maximum/minimum that are allowed for the sensitivity
    of each species (individually)    
    
    av means average, t means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data 6"""
    
    ave = uni(ave_min,ave_max, num) # average sensitivity
    alpha = uni(-0.95,-0.05, num) #interaction coefficient
    
    # does not affect coexistence
    t_f = uni(-1/sqrt, 1/sqrt,num) # stdv/mean of per capita contribution

    fix = np.array(num*[False])
    not_fix = np.logical_not(fix)
    t_e_fix = 100000*np.ones(num)
    t_mu_fix = 100000*np.ones(num)
    n_fix = 100000*np.ones(num)
    # runs until community that coexists is found
    # WARNING: changing settings might turn this into an infinite loop
    while num>0:
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        
        n = np.random.randint(5,21,num) # number of species
        
        comp = -alpha[not_fix]*n/(1-alpha[not_fix]*(n-1)) #effective competition
        t_mu = vec_uni(comp-1,1-comp)/sqrt
        minimum = np.amin(np.abs([e_max/ave[not_fix]-1, 1-e_min/ave[not_fix]]), axis = 0)
        #ensures, that sensitivity e_i  is in [e_min, e_max]
        t_e = vec_uni(-minimum/sqrt,minimum/sqrt)
        
        coex = coex_test(ave[not_fix],t_e,t_mu,comp)
        
        
        t_e_fix[not_fix] = t_e
        t_mu_fix[not_fix] = t_mu
        n_fix[not_fix] = n
        fix[not_fix] = coex
        not_fix = np.logical_not(fix)
        num = len(t_e_fix[not_fix])
    return ave,t_e_fix,t_mu_fix,\
            t_f,-alpha*n_fix/(1-alpha*(n_fix-1)),alpha,n_fix

def coex_test(ave,t_e,t_mu,comp):
    """tests if coexistence is given in changed site
    
    returns True, iff coexistence is guaranteed, see supplemental Info 11"""
    #computes the growthrate of species in the changed site
    mu_change = lambda u_i: (1+u_i*t_mu*sqrt)*(1/ave-(1+u_i*t_e*sqrt))
    #miminum always occurs on boundaries
    min1 = mu_change(1)/(1/ave-1-t_mu*t_e)
    min2 = mu_change(-1)/(1/ave-1-t_mu*t_e)
    return np.amin([min1, min2], axis = 0)>comp
    
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
