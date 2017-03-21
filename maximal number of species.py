import numpy
import matplotlib.pyplot as plt
"""computes and plots the number of species for a maximal EF,
depending on rho, the covariance and alpha, the interaction coefficient"""


""" plots the function between x_min, x_max"""
def plotter(x_min,x_max,function,accuracy=1000):
    x_values=numpy.linspace(x_min,x_max,accuracy)
    
    y_values=[function(x) for x in x_values]
    plt.plot(x_values,y_values,label='alpha='+str(n))
        

""" this computes the maximal n"""            
def function(rho):
    return   numpy.round((2+alpha-numpy.sqrt((2+alpha)**2-4*(1+alpha)*(1+(1+alpha)/rho)))/(2*alpha),0)
    
def function2(rho):
    return  numpy.round(1+(rho+1+alpha)/rho/alpha)
    


maximum=0
rho_max=0
x_max=0 
for alpha in numpy.linspace(-1,-0.2,10):
    plotter(-1,-0.15,function)
    #plotter(-1,-0.15,function2)
    for rho in numpy.linspace(-1,-0.2,10):
        save=function(rho)
        if save>maximum:
            maximum=save
            rho_max,x_max=rho,x
        
print("maximum is ",maximum)

x=x_max

