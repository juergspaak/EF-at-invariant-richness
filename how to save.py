"""Saving values"""
import numpy as np

"""save figures:
 fig = plt.figure()
plt.plot(range(10))
fig.savefig('temp.pdf', dpi=fig.dpi) #saves with same resolution
"""    

# Example
x = range(10)
y = range(10)
np.savez("savefile.npz", x = x,y = y)

#Can close now everything, next time:
    
npz = np.load("savefile.npz")
#shows all saved arrays
npz.files

npz05 = np.load("save_p05.npz")
npz25 = np.load("save_p25.npz")
npz50 = np.load("save_p50.npz")
npz75 = np.load("save_p75.npz")
npz95 = np.load("save_p95.npz")
npzrand = np.load("save_prand.npz")

np.savez("save_all.npz", paras_p05 = npz05['paras_p05'], para_p05 = npz05['para_p05'], 
         delta_EF_p05 = npz05['delta_EF_p05'], paras_p25 = npz25['paras_p25'], 
        para_p25 = npz25['para_p25'], delta_EF_p25 = npz25['delta_EF_p25'],
        paras_p50 = npz50['paras_p50'], para_p50 = npz50['para_p50'], 
        delta_EF_p50= npz50['delta_EF_p50'],paras_p75 = npz75['paras_p75'], 
        para_p75 = npz75['para_p75'], delta_EF_p75 = npz75['delta_EF_p75'],
        paras_p95 = npz95['paras_p95'], para_p95 = npz95['para_p95'], 
        delta_EF_p95 = npz95['delta_EF_p95'],paras_prand = npzrand['paras_prand'], 
        para_prand = npzrand['para_prand'], delta_EF_prand = npzrand['delta_EF_prand'],
        paras_p100 = paras_p100, para_p100 = para_p100, delta_EF_p100 = delta_EF_p100,
        keys05 = keys05, keys25 = keys25,keys50 = keys50,keys75 = keys75,
        keys95 = keys95,keysrand=keysrand, key_all = key_all, 
        key_relevant= key_relevant, keys_paras =keys_paras, keys_para= keys_para,
        keys_delta_EF =keys_delta_EF, labels = labels)

keys05 =['paras_p05','para_p05','delta_EF_p05']
keys25 =['paras_p25','para_p25','delta_EF_p25']
keys50 =['paras_p50','para_p50','delta_EF_p50']
keys75 =['paras_p75','para_p75','delta_EF_p75']
keys95 =['paras_p95','para_p95','delta_EF_p95']
keys100 =['paras_p100','para_p100','delta_EF_p100']
keysrand =['paras_prand','para_prand','delta_EF_prand']

keys_paras = ["paras_p05", "paras_p25", "paras_p50", "paras_p75", "paras_p95", "paras_prand"]
keys_para = ["para_p05", "para_p25", "para_p50", "para_p75", "para_p95", "para_prand"]
keys_delta_EF = ["delta_EF_p05", "delta_EF_p25", "delta_EF_p50", "delta_EF_p75", 
                 "delta_EF_p95", "delta_EF_prand","delta_EF_p100"]
                 
labels = ["p = 0.05","p = 0.25","p = 0.50","p = 0.75","p = 0.95","p random","p = 1.00"]
                 
key_all = keys05+keys25+keys50+keys75+keys95+keysrand
key_relevant = keys05+keys50+keys95+keys100



np.savez("save_p05.npz", paras_p05 = paras_p05, para_p05 = para_p05, delta_EF_p05 = delta_EF_p05)
np.savez("save_p25.npz", paras_p25 = paras_p25, para_p25 = para_p25, delta_EF_p25 = delta_EF_p25)
np.savez("save_p50.npz", paras_p50 = paras_p50, para_p50 = para_p50, delta_EF_p50 = delta_EF_p50)
np.savez("save_p75.npz", paras_p75 = paras_p75, para_p75 = para_p75, delta_EF_p75 = delta_EF_p75)
np.savez("save_p95.npz", paras_p95 = paras_p95, para_p95 = para_p95, delta_EF_p95 = delta_EF_p95)
np.savez("save_prand.npz", paras_prand = paras_prand, para_prand = para_prand, delta_EF_prand = delta_EF_prand)