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


sqrt = np.sqrt(3) #is needed often in the program
n = 20

def rand_par_repl(count=False,p = 'rand', ave_max = 0.5, e_max = 1,
                  ave_min = -0.5, e_min = -1):
    """ returns randomized parameters for one ecosystem
    av means average, t_ means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data, appendix G, folder results"""
    
    # generate distribution of per capita contributions for species types
    t_fc, t_fr, t_fs = uni(-1/sqrt, 1/sqrt,3) # stdv/mean
    avfc, avfr, avfs = uni(0.5,1.5,3) #averages of f
    

    counter = 0
    # fixed parameters, will only be changed if "no" community can be found
    # having these parameters
    avec = uni(ave_min,ave_max) # average effect on C species
    alpha = uni(-0.95,-0.05) # interaction coecfficient
    comp = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
    
    # percent of species with type C
    if p is 'rand':
        p = int(np.random.randint(1,n-1))/n
    q = 1-p
    
    #randomly generate communities, until one fullfills coexistence conditions
    #attention, changing settings might turn this into an infinite loop
    while not (counter and coex_test_repl(*parameters)):
        # chosen such that avmur,acmuc/(p*avmuc+q*avmuc)>comp
        avmuc = uni(0,10)
        avmur = uni(avmuc*comp*p/(1-q*comp),
                             min(avmuc*(1-p*comp)/(q*comp),10))
        
        
        treshhold1 = comp*(p*avmuc+q*avmur)
        # chosen such that min (mur,muc) > treshhold1, i.e. coexist
        t_muc = uni(-(1-treshhold1/avmuc)/sqrt,(1-treshhold1/avmuc)/sqrt)
        t_mur = uni(-(1-treshhold1/avmur)/sqrt,(1-treshhold1/avmur)/sqrt)
        
        # acmus*(1-ave_min)/(p*avmuc)(1-avec)+q*avmus(1-ave_min))>comp
        treshhold2 = avmuc*(1-avec)*p*comp/(1-comp*q)/(1-ave_min)
        # we always have treshhold2<treshhold1
        avmus = uni(treshhold2,treshhold1)
        
        # ensure, that min(mus) fullfills above conditions
        dist = min(treshhold1/avmus-1, 1-treshhold2/avmus)/sqrt
        t_mus = uni(-dist,dist)
        
        # choose min, s.t. avmus*(1-aves) fullfills coexistence conditions
        treshhold1 = max(1-avmuc/avmus*(1-avec)*(1-comp*p)/(q*comp),ave_min)
        # choose max, s.t. avmuc*(1-avec) fullfills coexistence conditions
        treshhold2 = min(1-avmuc/avmus*(1-avec)/(1-comp*q)*(p*comp),ave_max)
        aves = uni(treshhold1, treshhold2)
        
        # choose borders, that e_i are within [e_min, e_max]
        minimum = min(np.abs([e_max/avec-1, 1-e_min/avec]))
        t_ec = uni(-minimum/sqrt,minimum/sqrt)
        minimum = min(np.abs([e_max/aves-1, 1-e_min/aves]))
        t_es = uni(-minimum/sqrt,minimum/sqrt)
        
        # reference types are assumed to have e_i = 1, always         
        #aver = 1 #change if desired differently
        #t_er = 0
        if counter >10000: # avoids infinite loops, happens about 1/150
            counter = 0 
            avec = uni(0,0.5) # redefine new sensitivity
            alpha = uni(-0.95,-0.05) # interaction coecfficient
            comp =  -alpha*n/(1-alpha*(n-1))
        
        parameters = avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p
        counter+=1
    return  parameters # exited while loop, community fullfills coexistence
    
def coex_test_repl(avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p):
    """tests if coexistence is given in changed site
    
    input should be the output of rand_par_repl
    
    Note: Does not check coexistence conditions in reference site,
    nor that U species cannot survive in changed site and vice versa
    These conditions are always fullfilled by the chosen parameter settings
    By changing the above parameters this might become necessary"""
    
    
    q = 1-p
    sign_c, sign_s = np.sign([avec, aves]) #changes sign of inequalities
    
    # coex in stressed site
    avmuc_s = avmuc*avec*(1/avec-1 - t_muc*t_ec) #average growth rates of
    avmus_s = avmus*aves*(1/aves-1 - t_mus*t_es) #species in changed site
    #computes the growthrate of one species
    mu_str = lambda avmu, tmu, ave, te, u_i: \
                    avmu*(1+u_i*tmu*sqrt)*ave*(1/ave-(1+u_i*te*sqrt))
    # compute growthrates of extrema on boundaries:
    minimalc1 = mu_str(avmuc, t_muc, avec, t_ec,1) 
    minimalc2 = mu_str(avmuc, t_muc, avec, t_ec,-1)    
    minimals1 = mu_str(avmus, t_mus, aves, t_es,1)
    minimals2 = mu_str(avmus, t_mus, aves, t_es,-1)
    #minimal growthrate of all species in changed site
    minimal = min(minimalc1,minimalc2, minimals1, minimals2)\
                    /(p*avmuc_s+q*avmus_s)
                    
    """The following checks whether r species are extinct
    maximalr1 = mu_str(avmucr t_mur, aver, t_er,1) #maxima on the boundaries
    maximalr2 = mu_str(avmucr t_mur, aver, t_er,-1)
    # maxima in interior
    loc_r = 0.5*(1/t_er*(1/aver-1)-1/t_mur)
    if -1 < loc_r < 1:
        maximalr1 = max(maximalr1, 
            avmur*(1+loc_r*sqrt*t_mur)*aver*(1/aver-(1+loc_r*sqrt*t_er)))
    if max(maximalr1, maximalr2)/(p*avmuc_s+q*avmus_s)>comp
        return False"""
    
    if minimal<comp:
        return False
    else:
        return True

    
def EF(avmu1,avmu2,f1,f2,cov1,cov2,alpha,p, adjust):
    """ computes the EF of the given system
    
    For computational background see Eq. 6"""
    q = 1-p
    comp = -alpha*n/(1-alpha*(n-1))
    EF1 = n*f1*avmu1/(1+alpha)*(cov1+1-comp)
    EF2 = n*f2*avmu2/(1+alpha)*(cov2+1-comp)
    return p*EF1+q*EF2+adjust*p*q*n*comp/(1+alpha)*(f1-f2)*(avmu1-avmu2)
    
def rel_delta_EF_repl(avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p, adjust = 1):
    """computes \DeltaEF/EF in the case of changing composition
    
    For computational background see Eq. 7
    
    adjust can be set to 0 to see the effect of the adjustmentterms"""
    covc = t_muc*t_fc #covariances of the relative distributions
    covr = t_mur*t_fr
    EF_r = EF(avmuc,avmur,avfc,avfr,covc,covr,alpha,p, adjust)
    
    avmuc_s = avmuc*avec*(1/avec-1 - t_muc*t_ec) #average growth rates in
    avmus_s = avmus*aves*(1/aves-1 - t_mus*t_es) #changed site
    covc_s = t_fc*(t_muc*(1/avec-1)-t_ec)/(1/avec-1-t_ec*t_muc)
    covs_s = t_fs*(t_mus*(1/aves-1)-t_es)/(1/aves-1-t_es*t_mus)
    EF_s = EF(avmuc_s,avmus_s,avfc,avfs,covc_s,covs_s,alpha,p, adjust)

    return 100*(EF_s-EF_r)/EF_r #multiply by 100, result in percent

def rand_par_coex(ave_max = 0.5, e_max = 1,
                  ave_min = -0.5,e_min = -1):
    """ returns randomized parameters for one community"""
    
    ave = uni(ave_min,ave_max) # average sensitivity, in [0.01, 0.99]
    alpha = uni(-0.95,-0.05) #interaction coefficient
    parameters = [1,0,0,0,0,0,0] #values that do not coexist
    # runs until community that coexists is found
    #attention, changing settings might turn this into an infinite loop
    while not (coex_test_coex(*parameters)):
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        t_f = uni(-1/sqrt, 1/sqrt) # stdv/mean of per capita contribution, in [-1/sqrt(3),1/sqrt(3)]
        n = np.random.randint(5,21)
        comp = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
        t_mu = uni(comp-1,1-comp)/sqrt
        minimum = min(np.abs([e_max/ave-1, 1-e_min/ave]))
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

def asymptotic_function_repl(paras, max_ave_gam=1, test = False):
    if test:
        ave_gam = {'c':1000,'r':1000,'s':1000}
        t_gam = {'c':0,'r':0,'s':0}
    else:
        temp = uni(0,max_ave_gam,3)
        ave_gam = {'c':temp[0],'r':temp[1],'s':temp[2]}
        temp = uni(-1/sqrt, 1/sqrt,3)
        t_gam = {'c':temp[0],'r':temp[1],'s':temp[2]}
    #paras = rand_par_repl(ave_max = 0, e_max = 0)   
    avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p = paras
    avmu = {'c':avmuc,'r':avmur,'s':avmus}
    ave = {'c':avec,'r':1,'s':aves}
    avf = {'c':avfc,'r':avfr,'s':avfs}
    t_mu = {'c':t_muc,'r':t_mur,'s':t_mus}
    t_e = {'c':t_ec,'r':0,'s':t_es}
    t_f = {'c':t_fc,'r':t_fr,'s':t_fs}
    eco_fun = lambda x,t, N, H: avf[t]*(1+t_f[t]*x*sqrt)*H(x,t)*N(x,t)\
                                   /(N(x,t)+H(x,t))/2
    N = lambda x,t,mu,avmu: (mu(x,t)-comp*avmu)/(1+alpha)
    mu = lambda x,t: avmu[t]*(1+x*sqrt*t_mu[t])
    mu_change = lambda x,t: avmu[t]*(1+x*sqrt*t_mu[t])*\
                            (1-ave[t]*(1+t_e[t]*sqrt*x))
    N_ref = lambda x,t: N(x,t,mu,p*avmu['c']+(1-p)*avmu['r'])
    avmu_change = p*avmu['c']*(1-ave['c']*(1+t_mu['c']*t_e['c']))+\
                    (1-p)*avmu['s']*(1-ave['s']*(1+t_mu['s']*t_e['s']))
    N_change = lambda x,t: N(x,t,mu_change,avmu_change)
    H = lambda x,t: ave_gam[t]*(1+t_gam[t]*sqrt*x)*avmu[t]*(1+t_mu[t]*x*sqrt)
    EF = n*(p*quad(lambda x: eco_fun(x,'c',N_ref,H),-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'r',N_ref,H),-1,1)[0])
    EF_change = n*(p*quad(lambda x: eco_fun(x,'c',N_change,H),-1,1)[0]\
            +(1-p)*quad(lambda x: eco_fun(x,'s',N_change,H),-1,1)[0])
    delta_EF_lin = rel_delta_EF_repl(*paras)
    
    return 100*(EF_change-EF)/EF,delta_EF_lin,paras    