import numpy as np

def dimensionless_xi(D,rho_t,rho_pr, v, uesc, H, rho_0):
    """
        see equation 2 from Kurokawa et al (2018)
        
        D - diameter of impactor
        rho_t - density of mars
        rho_pr - density of impactor / projectile
        v - impact velocity
        uesc - escape velocity?
        H - atmospheric scale height
        rho_0 is atmospheric density at surface
        
    """
    
    xi = D**3*rho_t*rho_pr*(v**2-uesc**2) / (H**3*rho_0*uesc**2*(rho_t+rho_pr))
    
    return xi

def normalised_mass_lost_atmos(xi):
    """
        see equation 4, Kurokawa et al (2018)
        normalised atmosphere mass lost 
    """
    
    logXa = -6.375+5.239*np.log10(xi)-2.121*np.log10(xi)**2 + \
        0.397*np.log10(xi)**3-0.037*np.log10(xi)**4+0.0013*np.log10(xi)**5
    
    return 10**(logXa)

def fractional_change_atmos(xi,v,uesc):
    """
        see equation 4, Kurokawa et al (2018)
        normalised atmosphere mass lost - there seems to be a typo in the Kurokawa paper
        as it divides by the Matm, but this should be the mass of the projectile
        according to shuvalov equation 4
    """
    
    frac1 = 10**(-6.375+5.239*np.log10(xi)-2.121*np.log10(xi)**2 + \
        0.397*np.log10(xi)**3-0.037*np.log10(xi)**4+0.0013*np.log10(xi)**5) * \
            (v**2-uesc**2)/uesc**2
    return frac1

def change_atmos(xi,v,uesc,mass):
    """
        see equation 4, Kurokawa et al (2018)
        normalised atmosphere mass lost - there seems to be a typo in the Kurokawa paper
        as it divides by the Matm, but this should be the mass of the projectile
        according to shuvalov equation 4
    """
    
    frac1 = 10**(-6.375+5.239*np.log10(xi)-2.121*np.log10(xi)**2 + \
        0.397*np.log10(xi)**3-0.037*np.log10(xi)**4+0.0013*np.log10(xi)**5) * \
            (v**2-uesc**2)/uesc**2*mass
    return frac1

def normalised_projectile_mass(rho_t,rho_pr,v,uesc,xi):
    """
        see equation 5, Kurokawa et al. (2018)
        normalised projectile mass lost 
    """
#     Xpr = np.min(np.array( \
#         [0.035*rho_t/rho_pr*v/uesc*(np.log10(xi)-1), 0.07*rho_t/rho_pr*v/uesc, 1.]))
    Xpr = np.minimum(np.minimum(0.035*rho_t/rho_pr*v/uesc*(np.log10(xi)-1),0.07*rho_t/rho_pr*v/uesc), \
         np.ones(np.shape(xi)))
    return Xpr
    
def escape_velocity(Ggrav,Mmars,Rmars):
    uesc=np.sqrt(2.*Ggrav*Mmars/Rmars)
    return uesc
    
if __name__=="__main__":
    # test dimensionless_xi
    rho_t=4000.
    rho_pr = 3300.
    D=1000.
    v=15000.
    Ggrav=6.67e-11
    Mmars=6.39e23
    Rmars=3389.5e3
    uesc=escape_velocity(Ggrav,Mmars,Rmars)
    H=5000.
    rho_0=1.
    
    xi=dimensionless_xi(D,rho_t,rho_pr, v, uesc, H, rho_0)
    print('Dimensionless parameter ' + str(xi))
    
    # test normalised_mass_lost_atmos
    Xa = normalised_mass_lost_atmos(xi)
    print('Normalised mass lost (atmos) ' + str(Xa))
    
    # test normalised_projectile_mass
    Xpr = normalised_projectile_mass(rho_t, rho_pr, v, uesc, xi)
    print('Normalised mass lost (projectile)' + str(Xpr))
    
    # more tests - see Figure 3, Shuvalov (2009):
    import matplotlib.pyplot as plt
    xi=np.logspace(np.log10(2),np.log10(10000000),1000)
    Xa = normalised_mass_lost_atmos(xi)
    fractional_change_in_mass = Xa*(v**2-uesc**2)/uesc**2 # see equation 3
    
    
    plt.ion() 
    fig=plt.figure(figsize=(14,5))
    plt.subplot(131)
    plt.plot(xi,Xa)
    plt.xscale('log')  
    plt.yscale('log')  
    plt.show()
    plt.ylabel(r'$X_a$')
    plt.xlabel(r'$\xi$')
    plt.title('(a) Normalised mass lost')

    plt.subplot(132)
    plt.plot(xi,fractional_change_in_mass)
    plt.xscale('log')  
    plt.yscale('log')  
    plt.show()
    plt.ylabel(r'$\frac{\Delta m}{m}$')
    plt.xlabel(r'$\xi$')
    plt.title('(b) Fractional change in mass')
    
    
    # figure 4 b
    Xpr = normalised_projectile_mass(rho_t, rho_pr, v, uesc, xi)
    plt.subplot(133)
    plt.plot(xi,Xpr)
    plt.xscale('log')  
    plt.ylabel('$X_{pr}$')
    plt.xlabel(r'$\xi$')
    plt.title('(c) Normalised projectile mass')
    
    fig.tight_layout()
    
    
    
