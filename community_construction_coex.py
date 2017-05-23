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
from scipy.integrate import simps

sqrt = np.sqrt(3)
        
def vec_uni(low, high):
    """returns a random variable between low and high, uniform distribution
    
    low and high must be array-like with same shape.
    low is the lower bound of the uniform distribution, hight is the upperbound
    
    returns: An array of same shape as low and high. Each entry in this array
    is a randomvariable uniformly distributed in low, high.
    
    Examples:
    low = np.array([0,0,2])
    high = np.array([1,2,4])
    result = vec_uni(low, high)
    # result[0] is uniformly distributed in [0,1]
    # result[1] is uniformly distributed in [0,2]
    # result[2] is uniformly distributed in [2,4]"""
    high = np.array(high)
    low = np.array(low)
    return low+(high-low)*uni(size = low.shape)  #linear transformation

def rand_par(e_min=-1, ave_min = -0.5, ave_max = 0.5, e_max = 1,num = 100000):
    """returns randomized parameters for one community
    
    num : number of species to be generated
    ave_max and ave_min: maximum/minimum that are allowed for the
    average sensitivity of each species type
    e_max, e_min: maximum/minimum that are allowed for the sensitivity
    of each species (individually)    
    
    av means average, t means stdv/average,
    to understand the equations consult supplement data 6
    
    returns:
        ave: Average sensitivity for the communitie
        t_e: stdv/mean of sensitivity
        t_mu: stdv/mean of growthrates
        t_f: stdv/mean of per capita contribution to funciton
        comp: effective competition
        alpha: interaction coeficient
        n: number of species"""
    #check input correctness
    if not(e_min<=ave_min<=ave_max<=e_max):
        raise InputError("Please sort the input: e_min<=ave_min<=ave_max<=e_max")
    if e_max>1:
        raise InputError("e_max>1, effects above 1 are not allowed")
    #these values are fixed from the beginning, do not change while running
    ave = uni(ave_min,ave_max, num) # average sensitivity
    alpha = uni(-0.95,-0.05, num) #interaction coefficient

    # does not affect coexistence
    t_f = uni(-1/sqrt, 1/sqrt,num) # stdv/mean of per capita contribution
    
    # not_fix = False means values will not change in future runs
    not_fix = np.array(num*[True]) 
    #store values for return, these values will change
    t_e_fix = np.ones(num) #sdv/men of environmental change effect
    t_mu_fix = np.ones(num) #sdv/mean of growthrate
    n_fix = np.ones(num) #number of species
    # runs until all communities coexist
    # WARNING: changing settings might turn this into an infinite loop
    while num>0: #num is the number of communities that haven't been fiexed yet
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        n = np.random.randint(5,21,num) # number of species
        comp = -alpha[not_fix]*n/(1-alpha[not_fix]*(n-1)) #effective competition
        t_mu = vec_uni(comp-1,1-comp)/sqrt #min(mu)/mean(mu)>comp
        
        #ensures, that sensitivity e_i  is in [e_min, e_max]
        minimum = np.amin([np.sign(ave[not_fix])*(e_max/ave[not_fix]-1),
                    np.sign(ave[not_fix])*(1-e_min/ave[not_fix])], axis = 0)
        t_e = vec_uni(-minimum/sqrt,minimum/sqrt)
        
        #save ALL new computed values in not fixed communities
        t_e_fix[not_fix] = t_e
        t_mu_fix[not_fix] = t_mu
        n_fix[not_fix] = n
        #test communities coexistence requirements
        coex = coex_test(ave[not_fix],t_e,t_mu,comp)
        # fix communities that fullfill coexistence requirements
        not_fix[not_fix] = np.logical_not(coex)
        num = np.count_nonzero(not_fix) #number of not fixed communities
    return ave,t_e_fix,t_mu_fix,\
            t_f,-alpha*n_fix/(1-alpha*(n_fix-1)),alpha,n_fix

def coex_test(ave,t_e,t_mu,comp):
    """tests if coexistence is given in changed site
    
    returns True, iff coexistence is guaranteed, see supplemental Info 11"""
    #computes the growthrate of species in the changed site
    mu_change = lambda u_i: (1+u_i*t_mu*sqrt)*(1-ave*(1+u_i*t_e*sqrt))
    #miminum always occurs on boundaries
    min1 = mu_change(1)
    min2 = mu_change(-1)
    # min(mu(1-e))/mean(mu(1-e))>comp
    return np.amin([min1, min2], axis = 0)/(1-ave-ave*t_mu*t_e)>comp
    
def delta_EF_lin(ave,t_e,t_mu,t_f,comp,alpha = None,n=None):
    """ computes DetalEF/EF for one community
    
    ave,t_e,t_mu,t_f,comp,alpha,n: Return values of rand_par
    alpha and n are not needed, but passed s.t. one can run 
    delta_EF_lin(rand_par())
    
    for computational background see Eq. 4"""
    save_1 = -(1+t_mu*t_e)*ave
    save_2 = t_f*(t_mu+t_e)/(1+t_mu*t_e)-t_f*t_mu
    save_3 = t_mu*t_f+1-comp
    return 100*save_1*(1+save_2/save_3)
    
def EF_lin(ave,t_e,t_mu,t_f,comp,alpha = None,n=None):
    """ computes DetalEF/EF for one community
    
    ave,t_e,t_mu,t_f,comp,alpha,n: Return values of rand_par
    alpha and n are not needed, but passed s.t. one can run 
    delta_EF_lin(rand_par())
    
    for computational background see Eq. 4"""
    

def delta_EF_asym(ave,t_e,t_mu,t_f,comp,alpha,n,max_ave_H=1):
    """computes the EF with asymptotic f, f(N) = f_i*H_i*N_i/(N_i+H_i)
    
    ave,t_e,t_mu,t_f,comp,alpha,n should be the return values of rand_par
    max_ave_H is the maximum value for average value for  H.
    H_i is uniformly distributed in [0,2*ave_H]
    alpha and n are not needed, but passed s.t. one can run 
    delta_EF_lin(rand_par())"""
    num = len(ave) #number of communities
    # choose distribution of H: H ~u[0,2*ave]
    ave_H = uni(0,max_ave_H,num)
    t_H = uni(-1/sqrt, 1/sqrt,num) #stdv/mean of H
    H = lambda x: ave_H*(1+t_H*sqrt*x) #H_i for each species in a community
    
    #asymptotic EF in N, EF(N) = f_i*H_i*N_i/(N_i+H_i)
    eco_fun = lambda x, N: n*(1+t_f*x*sqrt)*H(x)*N(x)/(N(x)+H(x))
    
    # computes the equilibrium densities of species N, in changed and ref site
    N_ref = lambda x: (1+t_mu*sqrt*x-comp)/(1+alpha)
    N_change = lambda x: ((1+x*t_mu*sqrt)*(1-ave*(1+t_e*sqrt*x))-\
                comp*(1-ave*(1+t_mu*t_e)))/(1+alpha)
    
    # integrate over all species for EF
    x_simp = np.array(num*[np.linspace(-1,1,51)]) #x_axes
    y_ref = eco_fun(x_simp.T, N_ref).T #y_values in ref
    y_change = eco_fun(x_simp.T, N_change).T #y values in changed site
    EF_ref = simps(y_ref,x_simp) 
    EF_change = simps(y_change,x_simp)
    return 100*(EF_change-EF_ref)/EF_ref #multiply by 100 for percent
    
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        expr -- input expression in which the error occurred
        msg  -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = msg

para = {'e<0': rand_par(ave_max = 0, e_max = 0),
        'e>0': rand_par(ave_min = 0, e_min = 0)}
