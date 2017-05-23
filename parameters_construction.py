"""
@author: J.W. Spaak
This files generates the communities and saves them
"""
import pickle
from timeit import default_timer as timer
from warnings import warn

def para_return(fun):
    warn("The file 'repl, com_para.p' has not been found. It will be "+
                 "generated automatically and saved into the current folder."+
                 "This file will contain communitiy parameters that will be used"+
                 "for future use. This will take some minutes")
    
    
    num_com = int(1e5)
    #generate communities for different p
    keys = ['0.50','0.95' , 'rand', '0.05']
    count = 0
    para = {}
    start = timer()
    for key in keys:
        try:
            p= float(key)
        except ValueError:
            p = key
        para['e>0,'+key] = fun(num_com,
                            p=p,ave_min=0,e_min=0)
        now = timer()-start
        count+=1
        print("About {}/8 has been computed so far \n".format(count)+
               "time needed so far:{} seconds\nFinish ".format(now)+
               "within approximately {} seconds".format(now/count*(8-count)))
        
        para['e<0,'+key] = fun(num_com,
                            p=p,ave_max=0,e_max=0)
        now = timer()-start
        count+=1
        print("About {}/8 has been computed so far \n".format(count)+
               "time needed so far:{} seconds\nFinish ".format(now)+
               "within approximately {} seconds".format(now/count*(8-count)))
    
    #save data for future use  
    pickle.dump(para, open("repl, com_para.p","wb"))

    return para
