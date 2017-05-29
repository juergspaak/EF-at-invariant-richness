"""
@author: J.W. Spaak
Contains all functions, that are needed to generate communities as well as
compute the EF of those communities

Contains the same functions as community_construction_coex, but does so for
the replacing community structure
"""
import numpy as np
import pickle

from numpy.random import uniform as uni
from scipy.integrate import simps

sqrt = np.sqrt(3) #is needed often in the program
n = 20 #number of species, constant for all cases
      
def vec_uni(low, high):
    """returns a random variable between low and high, uniform distribution
    
    Parameters:
        low: array like
            lower boundaries for distributions
        high: array like
            higher boundaries for distributions, must have same shape as low
    
    Returns:
        rand_var: array like
            Same shape as low and high. Contains the random variables
    
    Examples:
    low = np.array([0,0,2])
    high = np.array([1,2,4])
    result = vec_uni(low, high)
    # result[0] is uniformly distributed in [0,1]
    # result[1] is uniformly distributed in [0,2]
    # result[2] is uniformly distributed in [2,4]"""
    high = np.array(high) #convert to arrays
    low = np.array(low)
    return low+(high-low)*uni(size = low.shape)  #linear transformation

def rand_par(e_min=-1, ave_min = -0.5, ave_max = 0.5, e_max = 1,p='rand',
              ad_com = 0.005,num = 100000):
    """ returns randomized parameters for num_com communities
    
    The function randomly generates num_com*(1+ad_com) communities until it
    finds num_com communities fullfilling all coexistence requirements. This
    method slightly shifts the distribution of alpha and e towards 0. The 
    error in the distribution is smaller than ad_com   
    Pleasse refer to supplementary data 7 to understand the code
    
    Input:
        num_com: scalar
            number of species to be generated
        p:  scalar or string
            Percent of species that are in ref and changed site (Type b species)
            p*n must be an integer or p ='rand'
            'rand will randomly generate values for p
            Note: p = 1 is NOT equivalent to coex.rand_par, because of the
            coexistence requirements.
        ave_max, ave_min: scalar<1, ave_min<=ave_max
            the maximum/minimum that are allowed for the
            average sensitivity of each species type
        e_max, e_min: scalar<1, e_min<=ave_min,ave_max<=e_max
            the maximum/minimum that are allowed for the sensitivity
            of each species (individually)
        ad_com: scalar
            Proportion of additionally computed communities.
    
    returns:
        mu: dict
            Contains all growth rates and some associated values
        e:  dict
            contains the sensitivity
        comp: array
            Relative competition
        f:  dict
            contains the per capita contributions
        p_ret:  array
            Percent of species that are in ref and changed site (Type b species)
            if p was a scalar, then p_ret = p*np.ones(num_com).
        alpha: array
            competition parameter        
        """
    #check input correctness
    if not(e_min<=ave_min<=ave_max<=e_max):
        raise InputError("Please sort the input: e_min<=ave_min<=ave_max<=e_max")
    if e_max>1: #growth rate of that species would be negative
        raise InputError("e_max>1, effects above 1 are not allowed")
    if not (p=='rand' or p*n==int(p*n)): #a species goes extinct or not
        raise InputError("p must either be 'rand' or p*n must be an integer")
        
    #save the original number of communites   
    num_com = num
    #number of communities to construct
    num = int(np.ceil(num_com*(1+ad_com)))
    
    # fixed parameters, do not change to find communities
    e_fix = {'avb':uni(ave_min,ave_max,num), #average effect on species b
            'avc': np.zeros(num), #will be filled with data while running
            'tc': np.zeros(num),
            'tb': np.zeros(num)} 
    mu_fix = {'avb':np.zeros(num),#will be filled with data while running
            'avc': np.zeros(num),
            'tc': np.zeros(num),
            'avu': np.zeros(num),
            'tu': np.zeros(num),
            'tb': np.zeros(num)}
    alpha = uni(-0.95,-0.05,num) # interaction coecfficient
    comp_fix = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
    # Fixed communities fullfill coexistence requirements
    not_fix = np.array(num*[True])
    # percent of species with type C
    if p is 'rand':
        p_fix = np.random.randint(1,n-1,num)/n
    else:
        p_fix = p*np.ones(num)
        
    #randomly generate communities, until num_com many fullfill coex. req.
    #attention, changing settings might turn this into an infinite loop
    while num>num_com*ad_com:
        #copy the predefined values into arrays to be used
        e = {'avb': e_fix['avb'][not_fix]}
        p = p_fix[not_fix]
        q = 1-p
        comp = comp_fix[not_fix]

        # min(mu['avb'],mu['avu'])/(p*mu['avb']+q*mu['avu'])>comp
        mu = {'avb': uni(0,10,num)}
        mu['avu'] = vec_uni(mu['avb']*comp*p/(1-q*comp),
            np.amin([mu['avb']*(1-p*comp)/(q*comp),10*np.ones(num)], axis = 0))
        
        #coexistence limit, min(mu)/mean(mu)>comp
        tresh1 = comp*(p*mu['avb']+q*mu['avu'])
        # chosen such that min (mu_u,mu_b) > tresh1, i.e. coexist
        mu['tb'] =vec_uni(-(1-tresh1/mu['avb'])/sqrt,(1-tresh1/mu['avb'])/sqrt)
        mu['tu'] =vec_uni(-(1-tresh1/mu['avu'])/sqrt,(1-tresh1/mu['avu'])/sqrt)
        
        # mu['avc']*(1-ave_min) must be able to coexist in changed site
        tresh2 = mu['avb']*(1-e['avb'])*p*comp/(1-comp*q)/(1-ave_min)
        # we always have treshhold2<treshhold1
        mu['avc'] = vec_uni(tresh2,tresh1)
        
        # ensure, that min(mu_c) fullfills same conditions as mu['avc']
        dist = np.amin([tresh1/mu['avc']-1, 1-tresh2/mu['avc']],axis = 0)/sqrt
        mu['tc'] = vec_uni(-dist,dist)
        
        # mu['avc']*(1-e['avc']) fullfills coexistence conditions
        # choose min for e['avc']
        tresh1 = np.amax([1-mu['avb']/mu['avc']*(1-e['avb'])*(1-comp*p)\
                     /(q*comp),ave_min*np.ones(num)],axis = 0)
        # choose max for e['avc']
        tresh2 = np.amin([1-mu['avb']/mu['avc']*(1-e['avb'])/(1-comp*q)\
                     *(p*comp),ave_max*np.ones(num)],axis = 0)
        e['avc'] = vec_uni(tresh1, tresh2)
        
        # choose borders, that e_i are within [e_min, e_max]
        minimum = np.amin([np.sign(e['avb'])*(e_max/e['avb']-1),
                np.sign(e['avb'])*(1-e_min/e['avb'])], axis = 0)
        e['tb'] = uni(-minimum/sqrt,minimum/sqrt)
        minimum = np.amin([np.sign(e['avc'])*(e_max/e['avc']-1),
                np.sign(e['avc'])*(1-e_min/e['avc'])], axis = 0)
        e['tc'] = vec_uni(-minimum/sqrt,minimum/sqrt)
        
        # average growthsrates in changed site of the species types
        mu['avb_change'] = mu['avb']*e['avb']*(1/e['avb']-1 - mu['tb']*e['tb'])
        mu['avc_change'] = mu['avc']*e['avc']*(1/e['avc']-1 - mu['tc']*e['tc'])
        # average growthrate of entire community in changed site
        mu['av_change'] = p*mu['avb_change']+q*mu['avc_change']
        
        # reference types are assumed to have e_i = 1, always 
        # if this part of the code is changed, please also change in coex_test        
        # e['avu'] = 1 #change if desired differently
        # e['tu'] = 0
        
        #copy the parameters into the fixed parameters
        for k in e_fix.keys():
            if k == 'avb': #do not copy into fixed 'avb'
                pass
            e_fix[k][not_fix] = e[k]
        for k in mu_fix.keys():
            mu_fix[k][not_fix] = mu[k]

        #check which species can coexist and update not_fix
        coex = coex_test(mu,e,comp)
        not_fix[not_fix] = np.logical_not(coex)
        num = np.count_nonzero(not_fix) #number of not fixed communities
        
    fix = np.logical_not(not_fix) #communities that are fixed, i.e. coex
    # choose only num_com coexisting communities
    comp_ret = comp_fix[fix][:num_com]
    alpha_ret = alpha[fix][:num_com]
    p_ret = p_fix[fix][:num_com]
    mu_ret = {key: mu_fix[key][fix][:num_com] for key in mu_fix.keys()}
    e_ret = {key: e_fix[key][fix][:num_com] for key in e_fix.keys()}
             
    # average growthsrates in changed site of the species types
    mu_ret['avb_change'] = mu_ret['avb']*e_ret['avb']*\
                        (1/e_ret['avb']-1 - mu_ret['tb']*e_ret['tb'])
    mu_ret['avc_change'] = mu_ret['avc']*e_ret['avc']*\
                        (1/e_ret['avc']-1 - mu_ret['tc']*e_ret['tc'])
    # average growthrate of entire community
    mu_ret['av_change'] = p_ret*mu_ret['avb_change']+\
                            (1-p_ret)*mu_ret['avc_change']
    
    # generate distribution of per capita contributions for species types
    t_fb, t_fu, t_fc = uni(-1/sqrt, 1/sqrt,[3,num_com]) # stdv/mean
    avfb, avfu, avfc = uni(0.5,1.5,[3,num_com]) #averages of f
    f = {'avb':avfb,'avu':avfu,'avc':avfc,\
         'tb':t_fb,'tu':t_fu,'tc':t_fc}
    # communities fullfill coexistence
    return mu_ret, e_ret,comp_ret,f,p_ret,alpha_ret
     

def coex_test(mu, e,comp, f = None, p = None, alpha = None):
    """tests if coexistence is given in changed site; see supp. Info 7
    
    Input: 
        mu, e, comp:
            As in output of rand_par
        f, alpha, p: optional
            Are not needed. They are just used s.t. one can
            run coex_test(rand_par(*args))
    
    returns:
        coex: array, dtype = boolean
            An array with coex rewuirements. True means, that this comunitiy
            fullfills the coexistence requirements
            
    Note: Does not check coexistence conditions in reference site,
    nor that U species cannot survive in changed site and vice versa
    These conditions are always fullfilled by the chosen parameter settings
    By changing the above parameters this might become necessary"""

    #computes the growthrate of one species
    mu_change = lambda x,t: mu['av'+t]*(1+x*mu['t'+t]*sqrt)\
                *e['av'+t]*(1/e['av'+t]-(1+x*e['t'+t]*sqrt))/mu['av_change']

    # minimal growthrate of all species in changed site, extremum is on boundary
    minimal = np.ones(len(comp))
    for x,t in [[1,'b'],[-1,'b'],[1,'c'],[-1,'c']]:
        minimal = np.amin([mu_change(x,t),minimal], axis = 0)
                
    """The following checks whether u species are extinct
    # maxima on the boundaries
    maximal = np.amax(mu_str(1,'u'),mu_str(-1,'u'),axis = 0) 
    # maxima in interior
    loc_max = 0.5*(1/e['tu']*(1/e['avu']-1)-1/mu['tu'])
    in_int = np.logical_and(-1<loc_max,loc_max<1)
    maximal[in_int]= np.amax(maximal[in_int],mu_str(loc_max[in_int],'u',axis=0)
    return np.logical(minimal>comp, maximal<comp)"""
    return minimal>comp

def EF_fun(mu,f,alpha,p,site,cov,adjust=True):
    """ computes the EF of the given system
    
    Input
        mu, f, alpha, p:
            as in output of rand_par
        s: "change" or "ref"
            Site EF is computed for
        cov: dict
            containing the covariances of mu and f
        adjust: boolean
            Set to False to see the effect f the adjusment terms       
    returns: 
        EF: array
            Array containing EF at site
        
    For computational background see Eq. 6"""
    s = {"ref": ['u',''], "change": ['c','_change']}[site]
    adjust = int(adjust)
    q = 1-p
    comp = -alpha*n/(1-alpha*(n-1))
    EF1 = n*f['avb']*mu['avb'+s[1]]/(1+alpha)*(cov['b'+s[1]]+1-comp)
    EF2 = n*f['av'+s[0]]*mu['av'+s[0]+s[1]]/(1+alpha)*(cov[s[0]+s[1]]+1-comp)
    return p*EF1+q*EF2+adjust*p*q*n*comp/(1+alpha)*(f['avb']-f['av'+s[0]])\
                            *(mu['avb'+s[1]]-mu['av'+s[0]+s[1]])

    
def delta_EF_lin(mu, e,comp,f,p, alpha,adjust = True):
    """computes \DeltaEF/EF in the case of changing composition
    
    For computational background see Eq. 7
    per capita contribution is assumed constant
    
    Input
        mu, e, comp, f, alpha, p:
            As in output of rand_par
        adjust: boolean
            Set to False to see the effect f the adjusment terms  
    returns: 
        deltaEF/EF: array
            Array containing 100*deltaEF/EF"""
    #covariances of the relative distributions
    cov = {'b': mu['tb']*f['tb'], 'u': mu['tu']*f['tu']} 
    cov['b_change'] = f['tb']*(mu['tb']*(1/e['avb']-1)-e['tb'])\
                        /(1/e['avb']-1-e['tb']*mu['tb'])
    cov['c_change'] = f['tc']*(mu['tc']*(1/e['avc']-1)-e['tc'])\
                        /(1/e['avc']-1-e['tc']*mu['tc'])
    #ecosystem funcitoning at reference site
    EF_u = EF_fun(mu,f,alpha,p,"ref",cov,adjust)
    #ecosystem funcitoning at changed site
    EF_c = EF_fun(mu,f,alpha,p,"change",cov,adjust)
    return 100*(EF_c-EF_u)/EF_u #multiply by 100, result in percent
    
def delta_EF_asym(mu, e,comp,f,p, alpha = None, max_ave_H=1):
    """computes the EF with asymptotic f, f(N) = f_i*H_i*N_i/(N_i+H_i)
    
    For more information see S10
    H_i is uniformly distributed in [0,2*ave_H]
    
    Input
        mu, e, comp, f, p:
            As in output of rand_par
        alpha: optional
            Is not needed. They are just used s.t. one can
            run delta_EF_asym
        max_ave_H: scalar, optional
            maximum for the average of H, maximum over all communities
        
    returns: 
        deltaEF/EF: array
            Array containing 100*deltaEF/EF, asymptotic contribution to EF"""
    num = len(alpha) #number of species
    # choose distributions of H: H ~u[0,2*ave]
    temp = uni(0,max_ave_H,3)
    gam = {'avb':temp[0],'avu':temp[1],'avc':temp[2]}
    temp = uni(-1/sqrt, 1/sqrt,3)
    gam.update({'tb':temp[0],'tu':temp[1],'tc':temp[2]})
    H = lambda x,t: gam['av'+t]*(1+gam['t'+t]*sqrt*x)\
                    *mu['av'+t]*(1+mu['t'+t]*x*sqrt)
    #asymptotic EF in N, f(N) = f_i*H_i*N_i/(N_i+H_i)
    eco_fun = lambda x,t, N: f['av'+t]*(1+f['t'+t]*x*sqrt)*H(x,t)*N(x,t)\
                                   /(N(x,t)+H(x,t))
    
    # growthrates in different sites
    mu_ref = lambda x,t: mu['av'+t]*(1+x*sqrt*mu['t'+t])
    mu_change = lambda x,t: mu['av'+t]*(1+x*sqrt*mu['t'+t])*\
                            (1-e['av'+t]*(1+e['t'+t]*sqrt*x))
    # computes the equilibrium densities of species N, in changed and ref site
    N = lambda x,t,mu,avmu: (mu(x,t)-comp*avmu)/(1+alpha)
    N_ref = lambda x,t: N(x,t,mu_ref,p*mu['avb']+(1-p)*mu['avu'])
    N_change = lambda x,t: N(x,t,mu_change,mu['av_change'])
    
    # integrate over all species for EF
    x_simp = np.array(num*[np.linspace(-1,1,51)]) #x_axes
    y_ref = {'b': eco_fun(x_simp.T, 'b',N_ref).T,
             'u': eco_fun(x_simp.T, 'u',N_ref).T}#y_values in ref
    y_cha = {'b': eco_fun(x_simp.T, 'b',N_change).T,
             'c': eco_fun(x_simp.T, 'c',N_change).T}#y_values in change
    # compute the EF
    EF_ref = n*(p*simps(y_ref['b'],x_simp)+(1-p)*simps(y_ref['u'],x_simp))
    EF_change = n*(p*simps(y_cha['b'],x_simp)+(1-p)*simps(y_cha['c'],x_simp))
    return 100*(EF_change-EF_ref)/EF_ref #multiply by 100 for percent
    
#load the communities
try:
    para = pickle.load(open("repl, com_para.p", "rb"))
except FileNotFoundError: #file not computed yet, will be computed
    import parameters_construction
    para = parameters_construction.para_return(rand_par)
    
class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class InputError(Error):
    """Exception raised for errors in the input.

    Attributes:
        msg  -- explanation of the error
    """

    def __init__(self, msg):
        self.msg = msg
