"""
@author: J.W. Spaak
Contains a non vectorized version of rand_par an coex_test

Contains the same functions as community_construction_coex, but does so for
the replacing community structure

"""

import numpy as np

from numpy.random import uniform as uni


sqrt = np.sqrt(3) #is needed often in the program
n = 20 #number of species in each community

def rand_par(p = 'rand', ave_max = 0.5, e_max = 1,
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
    while not (counter and coex_test(mu, e,f,comp,alpha,p)):
        if counter >10000: # avoids infinite loops, happens about 1/150
            counter = 1 
            e['avb'] = uni(ave_min,ave_max) # redefine new sensitivity
            alpha = uni(-0.95,-0.05) # interaction coecfficient
            comp =  -alpha*n/(1-alpha*(n-1))
            
        # chosen such that min(mu['avb'],mu['avu'])/(p*mu['avb']+q*mu['avu'])>comp
        mu = {'avb': uni(0,10)}
        mu['avu'] = uni(mu['avb']*comp*p/(1-q*comp),
                             min(mu['avb']*(1-p*comp)/(q*comp),10))
        
        #coexistence limit
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
        
        counter+=1 #count iterations to reset settings, avoids infinite loops
    return  mu, e,f,comp,alpha,p # community fullfills coexistence
     
  
def coex_test(mu, e,f,comp,alpha,p):
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
    minimal = 1
    for x,t in [[1,'b'],[-1,'b'],[1,'c'],[-1,'c']]:
        minimal = min(mu_change(x,t),minimal)
                
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
