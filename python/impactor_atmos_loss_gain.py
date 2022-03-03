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
    
    xi = D**3*rho_t*rho_pr*(v*2-uesc*2) / (H**3*rho_0*uesc**2*(rho_t+rho_pr))
    
    return xi

def normalised_mass_lost_atmos(xi):
    """
        see equation 4, Kurokawa et al (2018)
        normalised atmosphere mass lost 
    """
    
    logXa = -6.375+5.239*np.log(xi)-2.121*np.log(xi)**2 + \
        0.397*np.log(xi)**3-0.037*np.log(xi)**4+0.0013*np.log(xi)**5
    
    return np.exp(logXa)

def normalised_projectile_mass(rho_t,rho_pr,v,uesc,xi):
    """
        see equation 5, Kurokawa et al. (2018)
        normalised projectile mass lost 
    """
    
    Xpr = np.min(np.array( \
        [0.035*rho_t/rho_pr*v/uesc*(np.log(xi)-1), 0.07*rho_t/rho_pr*v/uesc, 1.]))
    return Xpr
    
    
if __name__=="__main__":
    # test dimensionless_xi
    rho_t=3000.
    rho_pr = 4500.
    D=1000.
    v=2000.
    uesc=100.
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
    
    
    