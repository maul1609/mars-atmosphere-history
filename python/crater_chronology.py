import numpy as np
# see Morbidelli et al. 2012 - assume projectile density of 2.6 g/m3 - see same reference

def N1(t):
    # see figure 1, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    return 5.44e-14*(np.exp(6.93*t)-1.)+8.38e-4*t
    
def dN1_dt(t):
    # see figure 1, Morbidelli et al. 2012
    # input time ago (i.e. should be positive)
    return 3.76992e-13*np.exp(6.93*t)+8.38e-4
    
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

def dN20_dt_integrand(t,model=1):
    if(model==1):
        dN=dN20_dt1(t)
    elif(model==2):
        if(t<4.1):
            dN=dN20_dt1(t)
        else:
            dN=dN20_dt2(t)
    elif(model==3):
        if(t<4.1):
            dN=dN20_dt1(t)
        else:
            dN=dN20_dt3(t)
    else:
        print('bad')
        dN=np.zeros(np.shape(t))
           
    return dN

def N1_between(Amars,t1,t2,model=1):
    from scipy.integrate import quad
    result = quad(dN1_dt,t1,t2)
    return Amars*result[0]

def N20_between(Amars,t1,t2,model=1):
    from scipy.integrate import quad
    result = quad(dN20_dt_integrand,t1,t2,model)
    return Amars*result[0]

if __name__ == "__main__":
    from scipy.integrate import quad
    
    Rmars=3389.5e3
    Amars=4.*np.pi*Rmars**2
    t=4.5
    # model 1 - just use  dN20_dt1 - blue line  
    result = quad(dN20_dt1,0,t)
    
    print('Model1: ' + str(Amars*result[0]))
    
    # model 2 - dN20_dt1 up to 4.1 and dN20_dt2 after - red line  
    result = quad(dN20_dt1,0,4.1)
    result2 = quad(dN20_dt2,4.1,t)
    
    print('Model2: ' + str(Amars*(result[0]+result2[0])))    
    
    # model 3 - dN20_dt1 up to 4.1 and dN20_dt3 after - red dashed line  
    result = quad(dN20_dt1,0,4.1)
    result2 = quad(dN20_dt3,4.1,t)
    
    print('Model3: ' + str(Amars*(result[0]+result2[0])))        
    
    
    # test overall
    result = quad(dN20_dt_integrand,0,t,2)
    print('Overall: ' + str(Amars*result[0]))
    