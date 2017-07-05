"""
@author: J.W. Spaak
This files generates the communities and saves them
"""
import pickle
from timeit import default_timer as timer
from warnings import warn

warn_message = "The file 'repl, com_para.p' has not been found. It will be "+\
    "generated automatically and saved into the current folder. "+\
    "This file will contain communitiy parameters to be used in the future. "+\
    "This will take some minutes."   
    
def para_return(fun):
    warn(warn_message)
    
    
    num = int(1e5)
    #generate communities for different, 
    keys = ['0.50','0.05','0.25', '0.75','0.95' ]
    count = 0
    para = {}
    start = timer()
    for key in keys:
        try:
            p= float(key)
        except ValueError:
            p = key
        para['e>0,'+key] = fun(p=p,ave_min=0,e_min=0, num = num)
        now = timer()-start
        count+=1
        print("About {}/10 have been computed so far \n".format(count)+
               "time needed so far:{} seconds\nFinish ".format(now)+
               "within approximately {} seconds\n".format(now/count*(10-count)))
        
        para['e<0,'+key] = fun(p=p,ave_max=0,e_max=0,num = num)
        now = timer()-start
        count+=1
        print("About {}/10 have been computed so far \n".format(count)+
               "time needed so far:{} seconds\nFinish ".format(now)+
               "within approximately {} seconds\n".format(now/count*(10-count)))
    
    #save data for future use  
    pickle.dump(para, open("repl, com_para.p","wb"))

    return para
