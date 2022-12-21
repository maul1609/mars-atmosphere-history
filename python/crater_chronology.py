import numpy as np
# see Morbidelli et al. 2012 - assume projectile density of 2.6 g/m3 - see same reference

def N1(t):
    # see figure 1, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    return 5.44e-14*(np.exp(6.93*t)-1.)+8.38e-4*t
    
def dN20_dt1(t):
    # see figure 3, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    N=5.9e-7+2.7e-16*np.exp(6.93*t)  

    return N

def dN20_dt2(t):
    # see figure 3, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    N=2.5e-2*np.exp(-((4.5-t)/0.003)**0.34)  

    return N

def dN20_dt3(t):
    # see figure 3, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    N=2.0e-2*np.exp(-((4.5-t)/0.01)**0.5)  

    return N


if __name__ == "__main__":
    from scipy.integrate import quad
    
    
    t=4.5
    # model 1 - just use  dN20_dt1 - blue line  
    result = quad(dN20_dt1,0,t)
    
    print('Model1: ' + str(result[0]))
    
    # model 2 - dN20_dt1 up to 3.5 and dN20_dt2 after - red line  
    result = quad(dN20_dt1,0,3.5)
    result2 = quad(dN20_dt2,3.5,t)
    
    print('Model2: ' + str(result[0]+result2[0]))    
    
    # model 3 - dN20_dt1 up to 3.5 and dN20_dt3 after - red dashed line  
    result = quad(dN20_dt1,0,3.5)
    result2 = quad(dN20_dt3,3.5,t)
    
    print('Model3: ' + str(result[0]+result2[0]))        