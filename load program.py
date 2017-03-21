"""
loads all variables"""

save_all = np.load("save_all.npz")
a = save_all

delta_EF_p05 = a['delta_EF_p05'].copy()
delta_EF_p25 = a['delta_EF_p25'].copy()
delta_EF_p50 = a['delta_EF_p50'].copy()
delta_EF_p75 = a['delta_EF_p75'].copy()
delta_EF_p95 = a['delta_EF_p95'].copy()
delta_EF_prand = a['delta_EF_prand'].copy()
print("one third")
para_p05 = a['para_p05'].copy()
para_p25 = a['para_p25'].copy()
para_p50 = a['para_p50'].copy()
para_p75 = a['para_p75'].copy()
para_p95 = a['para_p95'].copy()
para_prand = a['para_prand'].copy()
print("two thirds")
paras_p05 = a['paras_p05'].copy()
paras_p25 = a['paras_p25'].copy()
paras_p50 = a['paras_p50'].copy()
paras_p75 = a['paras_p75'].copy()
paras_p95 = a['paras_p95'].copy()
paras_prand = a['paras_prand'].copy()
print("done")
del a

a = np.load("all coex, 4e6 iterations.npz")
paras_p100 = a['paras'].copy()
delta_EF_p100 = a['rel_delta_EF'].copy()
del a