"""
@author: J.W. Spaak
This files generates the communities and saves them
"""
import pickle
from timeit import default_timer as timer

import help_functions as hef
import community_construction_repl as repl
import community_construction_coex as coex

#make an approximation, how long programm will run
start = timer()
hef.com_con(coex.rand_par,100,ave_min=0,e_min=0)
midpoint = timer()
hef.com_con(repl.rand_par,100,p='rand',ave_min=0,e_min=0)
hef.com_con(repl.rand_par,100,p='rand',ave_max=0,e_max=0)
end = timer()
num_com = int(1e5) #number of communities computed
approx_time_coex = num_com/100*2*(midpoint-start)
approx_time_repl = num_com/100*4*(end-midpoint)
print("This programm will take approximately {} seconds to terminate\n".format(
      approx_time_coex+approx_time_repl))

start = timer()
#generate communities, invariant species composition
com_para_coex = {}
com_para_coex['e>0'] = hef.com_con(coex.rand_par,num_com,ave_min=0,e_min=0)
com_para_coex['e<0'] = hef.com_con(coex.rand_par,num_com,ave_max=0,e_max=0)
#save data for future use
pickle.dump(com_para_coex, open("coex, com_para.p","wb"))

print("Progress report: Completed generation of invariant community structure."
      + "\nRuntime so far: "+str(timer()-start)
      + "\nApproximated time so far: " +str(approx_time_coex)
      + "\nApproximate time needed for the rest: " +str(approx_time_repl)+"\n")
#generate communities for different p
keys = ['0.95', '0.50', 'rand', '0.05']
count = 0
com_para_repl = {}
for key in keys:
    try:
        p= float(key)
    except ValueError:
        p = key
    com_para_repl['e>0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_min=0,e_min=0)
    #proportion that has already been done
    count+=1/8
    approx_time_sofar = approx_time_coex+count*approx_time_repl
    approx_time_rest = (1-count)*approx_time_repl
    print("Progress report: Completed {}, {},".format('e>0', key)
      + "\nRuntime so far: "+str(timer()-start)
      + "\nApproximated time so far: " +str(approx_time_sofar)
      + "\nApproximate time needed for the rest: " +str(approx_time_rest)+"\n")
    com_para_repl['e<0,'+key] = hef.com_con(repl.rand_par,num_com,
                        p=p,ave_max=0,e_max=0)
    count+=1/8
    approx_time_sofar = approx_time_coex+count*approx_time_repl
    approx_time_rest = (1-count)*approx_time_repl
    print("Progress report: Completed {}, {},".format('e<0', key)
      + "\nRuntime so far: "+str(timer()-start)
      + "\nApproximated time so far: " +str(approx_time_sofar)
      + "\nApproximate time needed for the rest: " +str(approx_time_rest)+"\n")

#save data for future use  
pickle.dump(com_para_repl, open("repl, com_para.p","wb"))
