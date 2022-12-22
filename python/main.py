import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import crater_chronology
import impactor_model
import impactor_velocity
import impactor_atmos_loss_gain_eq4_5

"""
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
"""
def test(y,t,omega):
    y1, z = y # this is y and dy/dx
    dydt = [z,-omega**2*y1 ]
    return dydt



if __name__ == "__main__":
    """
        initial conditions++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    Patm = 1.5 # initial pressure in most calculations - section 2.9, pp449
    tinit=-4.5e9 
    Rmars=3389.5e3 # m
    Mmars=6.39e23 # kg
    Ggrav=6.67e-11 # gravitational constant
    grav_mars=3.721 # m/s^2
    rho_t=4000. # density of martian surface
    tmars = 210. # temperature of martian surface - may change???
    Rgas = 8.314 # ideal gas constant
    MolW_atm = 44.01e-3 # molecular weight of martian atmosphere - may change???
    Hscale=Rgas*tmars/(MolW_atm*grav_mars) # scale height - m - may change???
    rho_pr=2600. # density of impactor
    
    total_impactor_mass = 2.e21 # kg
    tfinal=0.
    tstep=1e6
    omega=2.*np.pi*1.e-9
    output_interval = 1e7
    last_output = tinit-output_interval
    """
        ----------------------------------------------------------------------------------
    """

    """
        auxiliary calcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    Amars=4.*np.pi*Rmars**2 # surface area of mars
    Matm = Patm*1e5*Amars / grav_mars # mass of atmosphere
    t = np.mgrid[tinit:tfinal+tstep:tstep]
    nsteps = len(t)
    n_ode=3
    rand_numbers=np.random.rand((1000000))
    rand_numbers2=np.random.rand((1000000))
    """
        ----------------------------------------------------------------------------------
    """

    """
        set-up arrays ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    ystore=np.zeros((nsteps,n_ode))
    ystore[0,0]=1.
    ystore[0,1]=Patm
    # sample the sizes of the impactors
    diams, mass = impactor_model.sample_sizes(rand_numbers, rho_pr)
    # sample the impact velocities - not sure if this should be a different random number
    vels = impactor_velocity.inv_cdf1(rand_numbers2)
    # total number of impacts over the whole history
    num_impactor1_tot=crater_chronology.N1_between(Amars,0.,4.5)
    
    # escape velocity
    uesc=impactor_atmos_loss_gain_eq4_5.escape_velocity(Ggrav,Mmars,Rmars)
    
    """
        ----------------------------------------------------------------------------------
    """

    print(total_impactor_mass/np.sum(mass))
    #nsteps=0

    """
        main loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    for i in range(nsteps-1):
                
        
        
        # solve ODEs
        y0 = ystore[i,0:2]
        sol = odeint(test, y0, [t[i], t[i+1]], args=(omega,) )
        
        # 1. collapse? 0.5 bar and assumed to be 6 mbar
        
        # 2. calculate H2O - two temperatures
        
        # 3. impacts of asteroids and comets. during this time-step
        num_impactor1=crater_chronology.N1_between(Amars,-t[i+1]/1e9,-t[i]/1e9)
        num_impactor20=crater_chronology.N20_between(Amars,-t[i+1]/1e9,-t[i]/1e9,model=1)
        

        # note that results at each time get scaled by:
        # result * total_impactor_mass / np.sum(mass) * num_impactor1 / num_impactor1_tot
        scaling = total_impactor_mass / np.sum(mass) * num_impactor1 / num_impactor1_tot
        
        
        
        # 4. loss of the atmosphere following empirical formulae
        rho_0 = Patm*1.e5/tmars/(Rgas/MolW_atm)
        xi=np.maximum(impactor_atmos_loss_gain_eq4_5. \
            dimensionless_xi(diams,rho_t,rho_pr, vels, uesc, Hscale, rho_0)  , 1.e-10) 
        deltaM=impactor_atmos_loss_gain_eq4_5.change_atmos(xi,vels,uesc,mass)
        
        
        # left in the projectiles is (1-Xpr)*mass, but then only a certain % of this 
        #Xpr=np.maximum(impactor_atmos_loss_gain_eq4_5. \
        #    normalised_projectile_mass(rho_t,rho_pr,vels,uesc,xi),0)
        #deltaM=deltaM-Xpr*mass*0.01
        
        Matm = Matm - np.sum(deltaM)*scaling
        
        # 5. Sputtering and Photochemical escape
        
        # 6. Volcanic degassing - see rates
        #Matm = Matm + 0.5e27*MolW_atm*tstep*86400*365/6.02e23
        # 7. Isotope fractionation. 
        
        
        
        
        Patm = Matm*grav_mars/Amars/1.e5
        ystore[i+1,0]=sol[-1,0]
        ystore[i+1,1]=Patm
        
        if ((t[i]-last_output)>=output_interval):
            print('time: ' + str(t[i]/1e9) + '; num impacts 1: ' + str(num_impactor1) + \
                '; num impacts 20: ' + str(num_impactor20) + '; ratio: ' + \
                 str(num_impactor1/num_impactor20) + '; Patm: ' + str(Patm) )
            last_output = t[i]

    print('Final time: ' + str(t[-1]/1e9))
    """
        ----------------------------------------------------------------------------------
    """
    
    
    

    """
    b=0.25
    c=5.0
    y0 = [np.pi -0.1, 0.0]
    t=np.linspace(0,10,101)

    sol = odeint(pend, y0,t, args=(b,c))
    
    import matplotlib.pyplot as plt
    plt.ion()
    plt.plot(t, sol[:,0], 'b', label='theta(t)')
    plt.plot(t, sol[:,1], 'g', label='omega(t)')
    plt.legend(loc='best')
    plt.grid()
    plt.show()
    """


