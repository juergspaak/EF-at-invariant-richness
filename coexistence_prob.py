# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 08:46:39 2017

@author: Jurg
"""
import community_construction as community
from timeit import default_timer as timer
iterations = 500000
print("estimated time:", iterations*1.5/10000)

start = timer()

acc = 100
dict_aves = [[] for i in np.linspace(0,0.5,endpoint=False, num = acc)]

aves, counts = np.zeros(iterations), np.zeros(iterations)

for i in range(iterations):
    counts[i],aves[i] = community.rand_par_coex(True)
    dict_aves[int(aves[i]*acc/0.5)].append(counts[i])
    
averages = [np.average(i) for i in dict_aves]
medians = [np.median(i) for i in dict_aves]

lens = [len(i) for i in dict_aves]
end = timer()
print(end-start)            
