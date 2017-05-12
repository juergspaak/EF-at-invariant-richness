"""
@author: J.W. Spaak
This program was used to make all plots, that use asymptotic functioning
"""
import help_functions as hef
import community_construction_coex as coex

#generate communities
num_com = int(1e5) #number of communities computed
com_para = {}
com_para['e>0'] = hef.com_con(coex.rand_par,num_com,ave_min=0,e_min=0)
com_para['e<0'] = hef.com_con(coex.rand_par,num_com,ave_max=0,e_max=0)

print("com_con")

keys = ['100','2', '1','0.5']
# contains all EF datas
EF_data = {}
#compute asymptotic functioning for different averages of H
for key in keys:
    print(key)
    EF_data["e<0,"+key] = \
        hef.comp_EF(coex.delta_EF_asym,com_para["e<0"],max_ave_H=float(key))
    EF_data["e>0,"+key] = \
        hef.comp_EF(coex.delta_EF_asym,com_para["e>0"],max_ave_H=float(key))
#compare to linear functioning        
EF_data["e<0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,com_para["e<0"])
EF_data["e>0,ref"] = \
        hef.comp_EF(coex.delta_EF_lin,com_para["e>0"])
keys.append("ref")
ticks = [-60,-40,-20,0,20,40,60,80,100]
labels = [r'$\overline{H}=100$',r'$\overline{H}=2$',r'$\overline{H}=1$',
          r'$\overline{H}=0.5$',"Ref. Linear functioning"]
fig,ax = hef.percentiles(EF_data, keys, y_max = 100, ticks = ticks,labels = labels,
                     color = ["red", "green", "blue", "black", "orange"])