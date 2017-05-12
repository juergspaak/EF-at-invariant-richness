"""
@author: J.W. Spaak
This program was used to make Fig.2
"""
import help_functions as hef
import community_construction_coex as coex



#generate communities
num_com = int(1e4) #number of communities computed
com_para = {}
com_para['e>0,'] = hef.com_con(coex.rand_par,num_com,ave_min=0,e_min=0)
com_para['e<0,'] = hef.com_con(coex.rand_par,num_com,ave_max=0,e_max=0)

print("com_con")

#split into different covariances
EF_covs = {'e<0,11':[], 'e<0,1-1':[], 'e<0,-11':[], 'e<0,-1-1':[],
           'e>0,11':[], 'e>0,1-1':[], 'e>0,-11':[], 'e>0,-1-1':[],}
for key in ['e>0,', 'e<0,']:
    for para in com_para[key]:
            a = str(int(np.sign(-para[1]*para[2])))
            b = str(int(np.sign(para[1]*para[3])))
            EF_covs[key+a+b].append(coex.delta_EF_lin(*para))

        
# plot results
ticks = [-60,-40,-20,0,20,40,60,80,100]
keys = ['11','-11','-1-1','1-1']
labels = [r'cov$(|e|,\mu) = 1$, cov$(e,f)=1$',
              r'cov$(|e|,\mu) = -1$, cov$(e,f)=1$',
            r'cov$(|e|,\mu) = -1$, cov$(e,f)=-1$',
            r'cov$(|e|,\mu) = 1$, cov$(e,f)=-1$']
color = ['purple','red', 'green', 'blue']
fig,ax = hef.percentiles(EF_covs, keys, y_max = 100, ticks = ticks,
        color = color, labels = labels)
