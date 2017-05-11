The files community_construction_coex, community_construction_repl contain
all the functions needed to generate communities and calculate their DeltaEF/EF.

To understand the calculations being done please have a look at:
community_construction_coex: Supplemental Infomation, Section 6 and equation 4 from the paper
community_construction_repl: Supplemental Infomation, Section 7 and equations 6,7 from the paper

Examples:
The examples are given for the community_construction_coex functions, but exactly the same would also work
for the _repl case.

################################################
#generate communities
com_para = {}
# a gerneral community, facing any sort of environmental change
com_para['general'] = rand_par_coex() #a general community
#a community facing nutrient enrichment, i.e. e<0
com_para['enrich'] = rand_par_coex(ave_max=0, e_max =0)
#a community facing toxic stress, i.e. e>0
com_para['tox'] = rand_par_coex(ave_min=0, e_min =0)
# a community facing temperature changes, most species will suffer, but some might provit
com_para['temp'] = rand_par_coex(ave_min=-0.2, e_min =-0.3)

#compute deltaEF/EF, linear functioning
delta_EF = {key:rel_delta_EF_coex(*com_para[key]) \
		for key in com_para.keys()}

#compute deltaEF/EF, asymptotic functioning
asym_deltaEF = {key:asymptotic_function_coex(*com_para[key]) \
		for key in com_para.keys()}

#compute deltaEF/EF, asymptotic functioning, with low halfsaturating constant
asym_deltaEF_lowH = {key:asymptotic_function_coex(*com_para[key], max_ave_H = 0.2) \
		for key in com_para.keys()}
################################################
Note comparing these 3 dictionaries makes not that much sense, because the result is mainly
driven by a random process. Comparing these dictionaries makes only sense with a much larger
set of communities
