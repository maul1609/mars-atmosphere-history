import numpy as np
import matplotlib.pyplot as plt



# start of function
def impactor_sd(D):
    """ 
        Kurokawa et al. (2018)
         Small (D<300 km) impactors have q=3.5
         Large (D>300 km) impactors have q=2
    """
    
    dNdD = np.zeros((len(D)))
    q = np.zeros((len(D)))
    a = np.ones((len(D)))
    
    ind,=np.where(D<300e3)
    q[ind] = 3.5
    ind,=np.where(D>=300e3)
    q[ind] = 2.
    
    
    # Curve matching for the impactor size distribution
    # 300e3**-3.5 = a * 300e3**-2
    # a=(300e3**-3.5)/(300e3**-2)
    # a = 300e3**(-3.5+2=-1.5)
    a[ind] = 300.e3**(-1.5)
        
    
    dNdD = a* D**-q
    
    return (dNdD,a)
    # end of function

def impactor_csd(D,normalisingFactor,a):
    """
        integrate SD from dmin to D 
    """
    if(D<300e3):
        integral1=((D**-2.5)/-2.5 - (1e-6**-2.5)/-2.5) 
        integral1 /= normalisingFactor
    else:
        integral1=((300e3**-2.5)/-2.5 - (1e-6**-2.5)/-2.5)
        integral2=a*(1./300e3-1./D) 
        integral1 += integral2
        integral1 /= normalisingFactor
        
        
    return integral1
        





D=np.logspace(-6,6,10000)
(imp,a)=impactor_sd(D)

# evaluate the integral of impactor_sd between 1e-6 and infinity
integral1=(300e3**-2.5)/-2.5 - (1e-6**-2.5)/-2.5

# evaluate the integral of impactor_sd between 1e-6 and inf
integral2=a[-1]/300e3

# total area under curve
integralt=integral1+integral2
# scale original distribution by area under curve
imp = imp / integralt


"""
    sample cumulative distribution to generate sizes of impactors
"""
rand_numbers=np.random.rand((1000000))
diams = np.zeros((1000000))
for i in range(len(rand_numbers)):
    diams[i] = 

# plt.ion()
# plt.show()
# plt.plot(D,imp)
# plt.yscale('log')





