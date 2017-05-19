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
from scipy.integrate import quad

#load the communities
"""try:
    para = pickle.load(open("repl, com_para.p", "rb"))
except FileNotFoundError: #other functions should be still be available
    def error():
        raise FileNotFoundError("No such file or directory: 'repl, com_para.p'"
            "\nPlease run the file 'parameters of communities.py "+
            "to create this file. Read the readme for further instructions.")
"""        
def vec_uni(low, high):
    try:
        return low+(high-low)*uni(size = len(low))
    except TypeError:
        return uni(low, high)

sqrt = np.sqrt(3) #is needed often in the program
n = 20
def rand_par(num = 1000, p = 'rand', ave_max = 0.5, e_max = 1,
                  ave_min = -0.5, e_min = -1):
    """ returns randomized parameters for one community
    
    av means average, t means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data 7
    
    ave_max and ave_min are the maximum/minimum that are allowed for the
    average sensitivity of each species type
    e_max, e_min are the maximum/minimum that are allowed for the sensitivity
    of each species (individually)   
    
    returns:
        mu: Contains all growth rates and some associated values
        f: contains the per capita contributions
        e: contains the sensitivity
        comp: Relative competition
        p: Percent of species that are in ref and changed site (Type B species)
        """
    
    # generate distribution of per capita contributions for species types
    t_fb, t_fu, t_fc = uni(-1/sqrt, 1/sqrt,[3,num]) # stdv/mean
    avfb, avfu, avfc = uni(0.5,1.5,[3,num]) #averages of f
    f = {'avb':avfb,'avu':avfu,'avc':avfc,\
         'tb':t_fb,'tu':t_fu,'tc':t_fc}
    
    #avoids infinite loops
    counter = 0
    # fixed parameters, will only be changed if "no" community can be found
    # having these parameters
    e_fix = {'avb':uni(ave_min,ave_max,num),
            'avc': np.zeros(num),
            'tc': np.zeros(num),
            'tb': np.zeros(num)} # average effect on C species
    e = {'avb': e_fix['avb']}
    mu_fix = {'avb':np.zeros(num),
            'avc': np.zeros(num),
            'tc': np.zeros(num),
            'avu': np.zeros(num),
            'tu': np.zeros(num),
            'tb': np.zeros(num)} # average effect on C species
    alpha_fix = uni(-0.95,-0.05,num) # interaction coecfficient
    alpha = alpha_fix
    comp_fix = -alpha_fix*n/(1-alpha_fix*(n-1)) #effective competition , computed
    comp = comp_fix
    not_fix = np.array(num*[True])
    # percent of species with type C
    if p is 'rand':
        p_fix = np.random.randint(1,n-1,num)/n
    else:
        p_fix = p*np.ones(num)
    p = p_fix
    q = 1-p
    
    #randomly generate communities, until one fullfills coexistence conditions
    #attention, changing settings might turn this into an infinite loop
    while num>0:
        if counter >10000: # avoids infinite loops, happens about 1/150
            counter = 1 
            e['avb'] = uni(ave_min,ave_max,num) # redefine new sensitivity
            alpha_fix[not_fix] = uni(-0.95,-0.05,num) # interaction coecfficient
            comp_fix[not_fix] =  -alpha_fix[not_fix]*n/(1-alpha_fix[not_fix]*(n-1))
            alpha = alpha_fix[not_fix]
            comp = comp_fix[not_fix]
            
        # chosen such that min(mu['avb'],mu['avu'])/(p*mu['avb']+q*mu['avu'])>comp
        mu = {'avb': uni(0,10,num)}
        mu['avu'] = vec_uni(mu['avb']*comp*p/(1-q*comp),
            np.amin([mu['avb']*(1-p*comp)/(q*comp),10*np.ones(num)], axis = 0))
        
        #coexistence limit
        tresh1 = comp*(p*mu['avb']+q*mu['avu'])
        # chosen such that min (mu_u,mu_b) > tresh1, i.e. coexist
        mu['tb'] =vec_uni(-(1-tresh1/mu['avb'])/sqrt,(1-tresh1/mu['avb'])/sqrt)
        mu['tu'] =vec_uni(-(1-tresh1/mu['avu'])/sqrt,(1-tresh1/mu['avu'])/sqrt)
        
        # mu['avc']*(1-ave_min) must be able to coexist in changed site
        tresh2 = mu['avb']*(1-e['avb'])*p*comp/(1-comp*q)/(1-ave_min)
        # we always have treshhold2<treshhold1
        mu['avc'] = vec_uni(tresh2,tresh1)
        
        # ensure, that min(mu_c) fullfills above conditions
        dist = np.amin([tresh1/mu['avc']-1, 1-tresh2/mu['avc']],axis = 0)/sqrt
        mu['tc'] = vec_uni(-dist,dist)

        # choose min, s.t. mu['avc']*(1-e['avc']) fullfills coexistence conditions
        tresh1 = np.amax([1-mu['avb']/mu['avc']*(1-e['avb'])*(1-comp*p)\
                     /(q*comp),ave_min*np.ones(num)],axis = 0)
        # choose max, s.t. mu['avc']*(1-e['avc']) fullfills coexistence conditions
        tresh2 = np.amin([1-mu['avb']/mu['avc']*(1-e['avb'])/(1-comp*q)\
                     *(p*comp),ave_max*np.ones(num)],axis = 0)
        e['avc'] = vec_uni(tresh1, tresh2)
        
        # choose borders, that e_i are within [e_min, e_max]
        minimum = np.amin(np.abs([e_max/e['avc']-1, 1-e_min/e['avb']]),axis =0)
        e['tb'] = uni(-minimum/sqrt,minimum/sqrt)
        minimum = np.amin(np.abs([e_max/e['avc']-1, 1-e_min/e['avc']]),axis =0)
        e['tc'] = vec_uni(-minimum/sqrt,minimum/sqrt)
        
        # average growthsrates in changed site of the species types
        mu['avb_change'] = mu['avb']*e['avb']*(1/e['avb']-1 - mu['tb']*e['tb'])
        mu['avc_change'] = mu['avc']*e['avc']*(1/e['avc']-1 - mu['tc']*e['tc'])
        # average growthrate of entire community
        mu['av_change'] = p*mu['avb_change']+q*mu['avc_change']
        
        # reference types are assumed to have e_i = 1, always         
        # e['avu'] = 1 #change if desired differently
        # e['tu'] = 0
        for k in e_fix.keys():
            e_fix[k][not_fix] = e[k]
        for k in mu_fix.keys():
            mu_fix[k][not_fix] = mu[k]
        coex = coex_test(mu,e,comp)
        not_fix[not_fix] = np.logical_not(coex)
        alpha = alpha_fix[not_fix]
        comp = comp_fix[not_fix]
        p = p_fix[not_fix]
        q = 1-p
        e = {'avb': e_fix['avb'][not_fix]}
        num = len(p)
        counter+=1 #count iterations to reset settings, avoids infinite loops
    mu_fix['avb_change'] = mu_fix['avb']*e_fix['avb']*(1/e_fix['avb']-1 - mu_fix['tb']*e_fix['tb'])
    mu_fix['avc_change'] = mu_fix['avc']*e_fix['avc']*(1/e_fix['avc']-1 - mu_fix['tc']*e_fix['tc'])
        # average growthrate of entire community
    mu_fix['av_change'] = p_fix*mu_fix['avb_change']+(1-p_fix)*mu_fix['avc_change']
    return  mu_fix, e_fix,f,comp_fix,alpha_fix,p_fix # community fullfills coexistence
     

def coex_test(mu, e,comp):
    """tests if coexistence is given in changed site; see supp. Info 11
    
    input should be the output of rand_par_repl
    
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
    maximal = max(mu_str(1,'u'),mu_str(-1,'u')) #maxima on the boundaries
    # maxima in interior
    loc_r = 0.5*(1/e['tu']*(1/e['avu']-1)-1/mu['tu'])
    if -1 < loc_r < 1:
        maximal = max(maximal, mu_str(loc,'u')
    if maximal/mu['av_change']>comp
        return False"""
    return minimal>comp
start = timer()
a = rand_par(100000, ave_max = 0, e_max =0)  
print(timer()-start)
def EF_fun(mu,f,alpha,p,s,cov,adjust=1):
    """ computes the EF of the given system
    
    For computational background see Eq. 6"""
    q = 1-p
    comp = -alpha*n/(1-alpha*(n-1))
    EF1 = n*f['avb']*mu['avb'+s[1]]/(1+alpha)*(cov['b'+s[1]]+1-comp)
    EF2 = n*f['av'+s[0]]*mu['av'+s[0]+s[1]]/(1+alpha)*(cov[s[0]+s[1]]+1-comp)
    return p*EF1+q*EF2+adjust*p*q*n*comp/(1+alpha)*(f['avb']-f['av'+s[0]])\
                            *(mu['avb'+s[1]]-mu['av'+s[0]+s[1]])

    
def delta_EF_lin(mu, e,f,comp,alpha,p, adjust = 1):
    """computes \DeltaEF/EF in the case of changing composition
    
    For computational background see Eq. 7
    
    adjust can be set to 0 to see the effect of the adjustmentterms"""
    #covariances of the relative distributions
    cov = {'b': mu['tb']*f['tb'], 'u': mu['tu']*f['tu']} 
    cov['b_change'] = f['tb']*(mu['tb']*(1/e['avb']-1)-e['tb'])\
                        /(1/e['avb']-1-e['tb']*mu['tb'])
    cov['c_change'] = f['tc']*(mu['tc']*(1/e['avc']-1)-e['tc'])\
                        /(1/e['avc']-1-e['tc']*mu['tc'])

    EF_u = EF_fun(mu,f,alpha,p,['u',''],cov,adjust)

    EF_c = EF_fun(mu,f,alpha,p,['c','_change'],cov,adjust)
    return 100*(EF_c-EF_u)/EF_u #multiply by 100, result in percent
    
def delta_EF_asym(mu, e,f,comp,alpha,p, max_ave_H=1):
    """computes the EF with asymptotic f, f(N) = f_i*H_i*N_i/(N_i+H_i)
    
    mu, e,f,comp,alpha,p should be the return values of rand_par_repl
    max_ave_H is the maximum value for average value for  H.
    H_i is uniformly distributed in [0,2*ave_H*mu['av']]"""
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
    
    # Integrate for the average, divide by 2 for normalisation
    EF_ref = n*(p*quad(lambda x: eco_fun(x,'b',N_ref)/2,-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'u',N_ref)/2,-1,1)[0])
    EF_change = n*(p*quad(lambda x: eco_fun(x,'b',N_change)/2,-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'c',N_change)/2,-1,1)[0])
    return 100*(EF_change-EF_ref)/EF_ref #multiply by 100 for percent

