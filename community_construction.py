# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 09:03:21 2017

@author: Jurg
Contains all functions that are needed to randomly generate a community
and calculate its deltaEF
"""

import numpy as np
import numpy.random as rand


sqrt = np.sqrt(3)
n = 20

def rand_par_repl(count=False):
    """ returns randomized parameters for one ecosystem
    av means average, t_ means stdv/average,
    t values are always between -1/sqrt(3), 1/sqrt(3)
    to understand the equations consult supplement data, appendix G, folder results"""
    
    rands = np.random.rand(6)
    t_fc = 2*rands[0]/sqrt-1/sqrt # t of per capita contribution, in [-1/sqrt(3),1/sqrt(3)]
    t_fr = 2*rands[1]/sqrt-1/sqrt 
    t_fs = 2*rands[2]/sqrt-1/sqrt 
    avfc = rands[3]+0.5             # average of per capita contribution
    avfr = rands[4]+0.5
    avfs = rands[5]+0.5
    
    counter1,counter2=0,0
    avec = rand.uniform(0,0.5)
    alpha = rand.uniform(-0.95,-0.05) # interaction coecfficient
    comp = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
    p = int(rand.uniform(0,n-1)+1)/n
    p = 0.95
    q = 1-p
    
    parameters = 0
    while not (counter1 and coex_test_repl(*parameters)):
        rands = rand.rand(13)
        avmuc = rand.rand()*10
        # chosen such that avmur,acmuc/(p*avmuc+q*avmuc)>comp
        avmur = rand.uniform(avmuc*comp*p/(1-q*comp),min(avmuc*(1-p*comp)/(q*comp),10))

        treshhold = comp*(p*avmuc+q*avmur)
        # chosen such that min mur,muc > treshhold, i.e. coexist
        t_muc = rand.uniform(-(1-treshhold/avmuc)/sqrt,(1-treshhold/avmuc)/sqrt)
        t_mur = rand.uniform(-(1-treshhold/avmur)/sqrt,(1-treshhold/avmur)/sqrt)
        # chosen such that max(mus) < treshhold, i.e. can't coexist
        
        treshhold2 = avmuc*(1-avec)*p*comp/(1-comp*q)
        safty = 0.0*(treshhold -treshhold2)/2
        avmus = rand.uniform(treshhold2+safty,treshhold-safty)
        dist = min(treshhold/avmus-1, 1-treshhold2/avmus)/sqrt
        t_mus = rand.uniform(-dist,dist)
        
        treshhold = max(1-avmuc/avmus*(1-avec)*(1-comp*p)/(q*comp),0)
        treshhold2 = min(1-avmuc/avmus*(1-avec)/(1-comp*q)*(p*comp),0.5)
        aves = rand.uniform(treshhold, treshhold2)
        # 0<ec<1
        t_ec = 2*rands[9]/sqrt*min(1,1/avec-1)- 1/sqrt*min(1,1/avec-1)
        t_es = 2*rands[10]/sqrt*min(1,1/aves-1)- 1/sqrt*min(1,1/aves-1)
        
        # reference types can't coexist            
        #aver = 1#1-treshhold*comp/avmur*rands[11]
        #t_er = 2*rands[12]/sqrt*min(1,1/aver-1)- 1/sqrt*min(1,1/aver-1)
        if counter1 >10000: # avoids infinite loops, happens about 1/150
            counter1 = 0 
            counter2 += 1
            avec = rand.uniform(0,0.5)
            alpha = rand.uniform(-0.95,-0.05) # interaction coecfficient
            comp =  -alpha*n/(1-alpha*(n-1))
            
        parameters = avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p
        counter1+=1
        
            
    if count:
        return 10000*counter2+counter1 
    else:
        return  parameters
    
def coex_test_repl(avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p):
    """tests if coexistence is given in stessed site"""
    # coex in reference site
    q = 1-p
    
    # coex in stressed site
    avmuc_s = avmuc*avec*(1/avec-1 - t_muc*t_ec)
    avmus_s = avmus*aves*(1/aves-1 - t_mus*t_es)
    treshhold = comp*(p*avmuc_s+q*avmus_s)
    # extrema on boundaries:
    minimalc1 = avmuc*(1+sqrt*t_muc)*avec*(1/avec-(1+sqrt*t_ec))
    minimalc2 = avmuc*(1-sqrt*t_muc)*avec*(1/avec-(1-sqrt*t_ec))
    minimals1 = avmus*(1+sqrt*t_mus)*aves*(1/aves-(1+sqrt*t_es))
    minimals2 = avmus*(1-sqrt*t_mus)*aves*(1/aves-(1-sqrt*t_es))
    """The following checks whether r species are extinct, only necessary if aver != 0
    maximalr1 = avmur*(1+sqrt*t_mur)*aver*(1/aver-(1+sqrt*t_er))
    maximalr2 = avmur*(1-sqrt*t_mur)*aver*(1/aver-(1-sqrt*t_er))
    # maxima in interior
    loc_r = 0.5*(1/t_er*(1/aver-1)-1/t_mur)
    if -1 < loc_r < 1:
        maximalr1 = max(maximalr1, 
                        avmur*(1+loc_r*sqrt*t_mur)*aver*(1/aver-(1+loc_r*sqrt*t_er)))"""
    
    minimal = min(minimalc1,minimalc2, minimals1, minimals2)
    #maximal = max(maximalr1, maximalr2)
    if minimal<treshhold:# or maximal>treshhold:
        return False
    else:
        return True

    
def EF(avmu1,avmu2,f1,f2,cov1,cov2,alpha,p, adjust):
    """ computes the EF of the given system"""
    q = 1-p
    comp = -alpha*n/(1-alpha*(n-1))
    EF1 = n*f1*avmu1/(1+alpha)*(cov1+1-comp)
    EF2 = n*f2*avmu2/(1+alpha)*(cov2+1-comp)
    return p*EF1+q*EF2+adjust*p*q*n*comp/(1+alpha)*(f1-f2)*(avmu1-avmu2)
    
def rel_delta_EF_repl(avmuc,t_muc, avmur, t_mur, avmus, t_mus, avec, \
              t_ec, aves,t_es,comp,t_fc,t_fr,t_fs, \
              avfc,avfr,avfs,alpha,p, adjust = 1):
    covc = t_muc*t_fc
    covr = t_mur*t_fr
    EF_r = EF(avmuc,avmur,avfc,avfr,covc,covr,alpha,p, adjust)
    
    avmuc_s = avmuc*avec*(1/avec-1 - t_muc*t_ec)
    avmus_s = avmus*aves*(1/aves-1 - t_mus*t_es)
    covc_s = t_fc*(t_muc*(1/avec-1)-t_ec)/(1/avec-1-t_ec*t_muc)
    covs_s = t_fs*(t_mus*(1/aves-1)-t_es)/(1/aves-1-t_es*t_mus)
    EF_s = EF(avmuc_s,avmus_s,avfc,avfs,covc_s,covs_s,alpha,p, adjust)
    if EF_r <0 or EF_s<0:
        print ("weird, EF_r or EF_s<0")

    return 100*(EF_s-EF_r)/EF_r

def rand_par_coex(count = False, max_alpha = -0.05, 
                  ave_min = 0.0,e_min = 0.0):
    """ returns randomized parameters for one community"""
    
    ave = rand.uniform(ave_min,0.5) # average sensitivity, in [0.01, 0.99]
    alpha = rand.uniform(-0.95,max_alpha)
    counter = 0 # counts the number of tries to find a community
    parameters = 0
    while not (counter and coex_test_coex(*parameters)):
        """to understand why these parameters are chosen in this way, 
        see appendix G in supplement data, folder results"""
        t_f = rand.uniform(-1/sqrt, 1/sqrt) # stdv/mean of per capita contribution, in [-1/sqrt(3),1/sqrt(3)]
        n = int(rand.uniform(5,21))
        comp = -alpha*n/(1-alpha*(n-1)) #effective competition , computed
        t_mu = rand.uniform(comp-1,1-comp)/sqrt
        minimum = min(np.abs([1/ave-1, 1-e_min/ave]))
        t_e = rand.uniform(-minimum/sqrt,minimum/sqrt)
        parameters = ave,t_e,t_mu,t_f,comp,alpha,n
        counter +=1

    if count:
        return counter, ave
    else: 
        return ave,t_e,t_mu,t_f,comp,alpha,n  #alpha and n are just passed for convenience

def coex_test_coex(ave,t_e,t_mu,t_f,comp,alpha,n):
    """tests if coexistence is given in stessed site
    returns True, iff coexistence is guaranteed
    min_coex tests the species with minimal mu_i (assuming t_mu>0)"""
    sign = np.sign(ave)
    
    min_coex = (1+sqrt*t_mu)*(1/ave-1-sign*sqrt*t_e) \
                /(1/ave-1-t_mu*t_e)
    max_coex = (1-sqrt*t_mu)*(1/ave-1+sign*sqrt*t_e) \
                /(1/ave-1-t_mu*t_e)
    return min_coex>comp and max_coex>comp 
    
def rel_delta_EF_coex(ave,t_e,t_mu,t_f,comp,alpha,n, adjust = True):
    """ computes \DetalEF/EF for one ecosystem
    save_2, save_3 are used for the exact value,
    save_2_, save_3_ are used assuming cov(f, mu)=0"""
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
    
    