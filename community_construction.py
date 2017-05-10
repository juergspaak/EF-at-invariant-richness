# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 09:03:21 2017

@author: Jurg
Contains all functions that are needed to randomly generate a community
and calculate its deltaEF
"""

import numpy as np
from numpy.random import uniform as uni
from scipy.integrate import quad
import community_construction as test


sqrt = np.sqrt(3) #is needed often in the program
n = 20

def rand_par_repl(count=False,p = 'rand', ave_max = 0.5, e_max = 1,
                  ave_min = -0.5, e_min = -1):
    """ returns randomized parameters for one ecosystem
    av means average, t means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data 7"""
    
    # generate distribution of per capita contributions for species types
    t_fb, t_fu, t_fc = uni(-1/sqrt, 1/sqrt,3) # stdv/mean
    avfb, avfu, avfc = uni(0.5,1.5,3) #averages of f
    f = {'avb':avfb,'avu':avfu,'avc':avfc,\
         'tb':t_fb,'tu':t_fu,'tc':t_fc}
    
    #avoids infinite loops
    counter = 0
    # fixed parameters, will only be changed if "no" community can be found
    # having these parameters
    e = {'avb':uni(ave_min,ave_max)} # average effect on C species
    alpha = uni(-0.95,-0.05) # interaction coecfficient
    comp = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
    
    # percent of species with type C
    if p is 'rand':
        p = int(np.random.randint(1,n-1))/n
    q = 1-p
    
    #randomly generate communities, until one fullfills coexistence conditions
    #attention, changing settings might turn this into an infinite loop
    while not (counter and coex_test_repl2(mu, e,f,comp,alpha,p)):
    # chosen such that min(mu['avb'],mu['avu'])/(p*mu['avb']+q*mu['avu'])>comp
        mu = {'avb': uni(0,10)}
        mu['avu'] = uni(mu['avb']*comp*p/(1-q*comp),
                             min(mu['avb']*(1-p*comp)/(q*comp),10))
        
        
        tresh1 = comp*(p*mu['avb']+q*mu['avu'])
        # chosen such that min (mu_u,mu_b) > tresh1, i.e. coexist
        mu['tb'] = uni(-(1-tresh1/mu['avb'])/sqrt,(1-tresh1/mu['avb'])/sqrt)
        mu['tu'] = uni(-(1-tresh1/mu['avu'])/sqrt,(1-tresh1/mu['avu'])/sqrt)
        
        
        # mu['avc']*(1-ave_min) must be able to coexist in changed site
        tresh2 = mu['avb']*(1-e['avb'])*p*comp/(1-comp*q)/(1-ave_min)
        # we always have treshhold2<treshhold1
        mu['avc'] = uni(tresh2,tresh1)
        
        # ensure, that min(mu_c) fullfills above conditions
        dist = min(tresh1/mu['avc']-1, 1-tresh2/mu['avc'])/sqrt
        mu['tc'] = uni(-dist,dist)
        
        # choose min, s.t. mu['avc']*(1-e['avc']) fullfills coexistence conditions
        tresh1 = max(1-mu['avb']/mu['avc']*(1-e['avb'])*(1-comp*p)\
                     /(q*comp),ave_min)
        # choose max, s.t. mu['avc']*(1-e['avc']) fullfills coexistence conditions
        tresh2 = min(1-mu['avb']/mu['avc']*(1-e['avb'])/(1-comp*q)\
                     *(p*comp),ave_max)
        e['avc'] = uni(tresh1, tresh2)
        
        # choose borders, that e_i are within [e_min, e_max]
        minimum = min(np.abs([e_max/e['avc']-1, 1-e_min/e['avb']]))
        e['tb'] = uni(-minimum/sqrt,minimum/sqrt)
        minimum = min(np.abs([e_max/e['avc']-1, 1-e_min/e['avc']]))
        e['tc'] = uni(-minimum/sqrt,minimum/sqrt)
        
        # average growthsrates in changed site of the species types
        mu['avb_change'] = mu['avb']*e['avb']*(1/e['avb']-1 - mu['tb']*e['tb'])
        mu['avc_change'] = mu['avc']*e['avc']*(1/e['avc']-1 - mu['tc']*e['tc'])
        # average growthrate of entire community
        mu['av_change'] = p*mu['avb_change']+q*mu['avc_change']
        
        # reference types are assumed to have e_i = 1, always         
        # e['avu'] = 1 #change if desired differently
        # e['tu'] = 0
        if counter >10000: # avoids infinite loops, happens about 1/150
            counter = 0 
            e['avb'] = uni(0,0.5) # redefine new sensitivity
            alpha = uni(-0.95,-0.05) # interaction coecfficient
            comp =  -alpha*n/(1-alpha*(n-1))
        counter+=1

    return  mu, e,f,comp,alpha,p # exited while loop, community fullfills coexistence
     
  
def coex_test_repl(mu, e,f,comp,alpha,p):
    """tests if coexistence is given in changed site
    
    input should be the output of rand_par_repl
    
    Note: Does not check coexistence conditions in reference site,
    nor that U species cannot survive in changed site and vice versa
    These conditions are always fullfilled by the chosen parameter settings
    By changing the above parameters this might become necessary"""

    #computes the growthrate of one species
    mu_str = lambda x,t: mu['av'+t]*(1+x*mu['t'+t]*sqrt)\
                *e['av'+t]*(1/e['av'+t]-(1+x*e['t'+t]*sqrt))/mu['av_change']

    # minimal growthrate of all species in changed site, extremum is on boundary
    minimal = 1
    for x,t in [[1,'b'],[-1,'b'],[1,'c'],[-1,'c']]:
        minimal = min(mu_str(x,t),minimal)
                
    """The following checks whether u species are extinct
    maximal = max(mu_str(1,'u'),mu_str(-1,'u')) #maxima on the boundaries
    # maxima in interior
    loc_r = 0.5*(1/e['tu']*(1/e['avu']-1)-1/mu['tu'])
    if -1 < loc_r < 1:
        maximal = max(maximal, mu_str(loc,'u')
    if maximal/mu['av_change']>comp
        return False"""

    if minimal<comp:
        return False
    else:
        return True
    
def EF_(mu,f,alpha,p,s,cov,adjust=1):
    """ computes the EF of the given system
    
    For computational background see Eq. 6"""
    q = 1-p
    comp = -alpha*n/(1-alpha*(n-1))
    EF1 = n*f['avb']*mu['avb'+s[1]]/(1+alpha)*(cov['b'+s[1]]+1-comp)
    EF2 = n*f['av'+s[0]]*mu['av'+s[0]+s[1]]/(1+alpha)*(cov[s[0]+s[1]]+1-comp)
    return p*EF1+q*EF2+adjust*p*q*n*comp/(1+alpha)*(f['avb']-f['av'+s[0]])\
                            *(mu['avb'+s[1]]-mu['av'+s[0]+s[1]])
    
def rel_delta_EF_repl(mu, e,f,comp,alpha,p, adjust = 1):
    """computes \DeltaEF/EF in the case of changing composition
    
    For computational background see Eq. 7
    
    adjust can be set to 0 to see the effect of the adjustmentterms"""
    cov = {'b': mu['tb']*f['tb'], 'u': mu['tu']*f['tu']} #covariances of the relative distributions
    cov['b_change'] = f['tb']*(mu['tb']*(1/e['avb']-1)-e['tb'])\
                        /(1/e['avb']-1-e['tb']*mu['tb'])
    cov['c_change'] = f['tc']*(mu['tc']*(1/e['avc']-1)-e['tc'])\
                        /(1/e['avc']-1-e['tc']*mu['tc'])

    EF_r = EF_2(mu,f,alpha,p,['u',''],cov,adjust)
    EF_s = EF_2(mu,f,alpha,p,['c','_change'],cov,adjust)
    return 100*(EF_s-EF_r)/EF_r #multiply by 100, result in percent

def rand_par_coex(ave_max = 0.5, e_max = 1,
                  ave_min = -0.5,e_min = -1):
    """ returns randomized parameters for one community"""
    
    ave = uni(ave_min,ave_max) # average sensitivity
    alpha = uni(-0.95,-0.05) #interaction coefficient
    parameters = [1,0,0,0,0,0,0] #values that do not coexist
    # runs until community that coexists is found
    # WARNING: changing settings might turn this into an infinite loop
    while not (coex_test_coex(*parameters)):
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        t_f = uni(-1/sqrt, 1/sqrt) # stdv/mean of per capita contribution, in [-1/sqrt(3),1/sqrt(3)]
        n = np.random.randint(5,21) # number of species
        comp = -alpha*n/(1-alpha*(n-1)) #effective competition
        t_mu = uni(comp-1,1-comp)/sqrt
        minimum = min(np.abs([e_max/ave-1, 1-e_min/ave]))
        #ensures, that sensitivity e_i  is in [e_min, e_max]
        t_e = uni(-minimum/sqrt,minimum/sqrt)
        parameters = ave,t_e,t_mu,t_f,comp,alpha,n

    return ave,t_e,t_mu,t_f,comp,alpha,n  #alpha and n are just passed for convenience

def coex_test_coex(ave,t_e,t_mu,t_f,comp,alpha,n):
    """tests if coexistence is given in changed site
    
    returns True, iff coexistence is guaranteed"""
    #computes the growthrate of species in the changed site    
    mu_str = lambda tmu, ave, te, u_i: \
                    (1+u_i*tmu*sqrt)*(1/ave-(1+u_i*te*sqrt))
    min1 = mu_str(t_mu, ave, t_e,1)/(1/ave-1-t_mu*t_e)
    min2 = mu_str(t_mu, ave, t_e,-1)/(1/ave-1-t_mu*t_e)
    return min(min1, min2)>comp
    
def rel_delta_EF_coex(ave,t_e,t_mu,t_f,comp,alpha,n):
    """ computes \DetalEF/EF for one ecosystem
    
    for computational background see Eq. 4"""
    save_1 = -(1+t_mu*t_e)*ave
    save_2 = t_f*(t_mu+t_e)/(1+t_mu*t_e)-t_f*t_mu
    save_3 = t_mu*t_f+1-comp
    return 100*save_1*(1+save_2/save_3)
    
def fun_covs(para, delta_EF):
    covs = {'11':[], '1-1':[], '-11':[], '-1-1':[]}
    for i in range(len(delta_EF)):
        a = str(int(np.sign(para[i][1]*para[i][2])))
        b = str(int(np.sign(para[i][1]*para[i][3])))
        covs[a+b].append(delta_EF[i])
    covs2 = {}
    covs2[r'cov$(e,\mu) = 1$, cov$(e,f)=1$'] = covs['11']
    covs2[r'cov$(e,\mu) = -1$, cov$(e,f)=-1$'] = covs['-1-1']
    covs2[r'cov$(e,\mu) = -1$, cov$(e,f)=1$'] = covs['-11']
    covs2[r'cov$(e,\mu) = 1$, cov$(e,f)=-1$'] = covs['1-1']
    return covs2
    
def asymptotic_function_coex(max_ave_gam=1, test = False):
    if test:
        ave_gam = 100
        t_gam = 0
    else:
        ave_gam = uni(0,max_ave_gam)
        t_gam = uni(-1/sqrt, 1/sqrt)
    ave,t_e,t_mu,t_f,comp,alpha,n = rand_par_coex(ave_max = 0, e_max = 0)
    int_fun = lambda x, N, H: n*(1+t_f*x*sqrt)*H(x)*N(x)/(N(x)+H(x))/2
    N = lambda x: (1+t_mu*sqrt*x-comp)/(1+alpha)
    H = lambda x: ave_gam*(1+t_gam*sqrt*x)*(1+t_mu*x*sqrt)
    N_change = lambda x: ((1+x*t_mu*sqrt)*(1-ave*(1+t_e*sqrt*x))-\
                comp*(1-ave*(1+t_mu*t_e)))/(1+alpha)
    EF = quad(lambda x: int_fun(x,N,H),-1,1)[0]
    EF_stress = quad(lambda x: int_fun(x,N_change,H),-1,1)[0]
    paras = [ave,t_e,t_mu,t_f,comp,alpha,n,ave_gam, t_gam]
    delta_EF_lin = rel_delta_EF_coex(*paras[:7])
    
    return 100*(EF_stress-EF)/EF,delta_EF_lin,paras

def asymptotic_function_repl(mu, e,f,comp,alpha,p
                              , max_ave_gam=1, test = False):
    if test:
        gam = {'avb':1000,'avu':1000,'avc':1000,\
               'tb':0,'tu':0,'tc':0}
    else:
        temp = uni(0,max_ave_gam,3)
        gam = {'avb':temp[0],'avu':temp[1],'avc':temp[2]}
        temp = uni(-1/sqrt, 1/sqrt,3)
        gam.update({'tb':temp[0],'tu':temp[1],'tc':temp[2]})

    eco_fun = lambda x,t, N, H: f['av'+t]*(1+f['t'+t]*x*sqrt)*H(x,t)*N(x,t)\
                                   /(N(x,t)+H(x,t))/2
    N = lambda x,t,mu,avmu: (mu(x,t)-comp*avmu)/(1+alpha)
    mu_ref = lambda x,t: mu['av'+t]*(1+x*sqrt*mu['t'+t])
    mu_change = lambda x,t: mu['av'+t]*(1+x*sqrt*mu['t'+t])*\
                            (1-e['av'+t]*(1+e['t'+t]*sqrt*x))
    N_ref = lambda x,t: N(x,t,mu_ref,p*mu['avb']+(1-p)*mu['avu'])
    N_change = lambda x,t: N(x,t,mu_change,mu['av_change'])
    avmu_change = p*mu['avb']*(1-e['avb']*(1+mu['tb']*e['tc']))+\
                    (1-p)*mu['avc']*(1-e['avc']*(1+mu['tc']*e['tc']))
    H = lambda x,t: gam['av'+t]*(1+gam['t'+t]*sqrt*x)\
                    *mu['av'+t]*(1+mu['t'+t]*x*sqrt)
    EF = n*(p*quad(lambda x: eco_fun(x,'b',N_ref,H),-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'u',N_ref,H),-1,1)[0])
    EF_change = n*(p*quad(lambda x: eco_fun(x,'b',N_change,H),-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'c',N_change,H),-1,1)[0])
    delta_EF_lin = rel_delta_EF_repl2(mu, e,f,comp,alpha,p)
    
    return 100*(EF_change-EF)/EF,delta_EF_lin
