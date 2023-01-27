import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import crater_chronology
import impactor_model
import impactor_velocity
import impactor_atmos_loss_gain_eq4_5
import svp_ice
import volcano_outgassing
import sputtering_co2

# table 4, Kurokawa et al.
isotopes = ['d15N',
            '20Ne/22Ne',
            '36Ar/38Ar',
            '78Kr/84Kr',
            '80Kr/84Kr',
            '82Kr/84Kr',
            '83Kr/84Kr',
            '86Kr/84Kr',
            '124Xe/130Xe',
            '126Xe/130Xe',
            '128Xe/130Xe',
            '129Xe/130Xe',
            '131Xe/130Xe',
            '132Xe/130Xe',
            '134Xe/130Xe',
            '136Xe/130Xe']
isotope_comps = dict()
isotope_comps['Mars'] = [572,10.1,4.2,6.37e-3,4.09e-2,0.2054,\
                         0.2034,0.3006,2.45e-2,2.12e-2,0.4767, \
                         16.400,5.147,6.460,2.587,2.294]
                         
isotope_comps['Volcanic degassing'] = [-30,13.7,5.3,5.962e-3,3.919e-2,\
                    0.20149,0.20141,0.30950,2.851e-2,2.512e-2,0.5073,6.358,\
                    5.043,6.150,2.359,1.988]
                    
isotope_comps['Asteroids'] = [-30, 8.9,5.3,5.962e-3,3.919e-2,\
                    0.20149,0.20141,0.30950,2.851e-2,2.512e-2,0.5073,6.358,\
                    5.043,6.150,2.359,1.988]
                    
isotope_comps['Comets'] = [1000,13.7,5.4,6.470e-3,4.124e-2,0.20629,0.20340,\
                        0.29915,2.947e-2,2.541e-2,0.50873,6.287,4.9958, \
                        6.0479,2.1288,1.6634]
                    
isotope_comps['IDPs'] = [np.nan, 13.7,5.8,6.470e-3,4.124e-2,0.20629,0.20340,\
                        0.29915,2.947e-2,2.541e-2,0.50873,6.287,4.9958, \
                        6.0479,2.1288,1.6634]


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
    #MolW_atm = 44.01e-3 # molecular weight of martian atmosphere - may change???
    MolW_atm = np.array([44.01e-3, 28.0134e-3, 18.01528e-3])
    rho_pr=2600. # density of impactor
    Pcollapse = 0.5
    sampleFlag = 0 # 0=once, 1=every time-step; 
    
    total_impactor_mass = 2.e21 # kg
    tfinal=0.
    tstep=1e6
    omega=2.*np.pi*1.e-9
    output_interval = 1e7
    last_output = tinit-output_interval
    
    isotope_comps_sim = np.array(isotope_comps['Volcanic degassing'])
    ind,=np.where(['Xe' in iso for iso in isotopes])
    isotope_comps_sim[ind]=np.array(isotope_comps['Mars'])[ind]
    """
        ----------------------------------------------------------------------------------
    """




    """
        auxiliary calcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    Amars=4.*np.pi*Rmars**2 # surface area of mars
    Matm = Patm*1e5*Amars / grav_mars # mass of atmosphere
    Matm = np.array([Matm, 0., 0.]) # mass of co2, n2, h2o
    Natm = Matm / MolW_atm     # number of moles of co2, n2, h2o
    ph2o = svp_ice.svp_ice01(tmars)
    Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
    t = np.mgrid[tinit:tfinal+tstep:tstep]

    nsteps = len(t)
    n_ode=3
    rand_numbers=np.random.rand((100000))
    rand_numbers2=np.random.rand((100000))
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
        if sampleFlag == 1:
            rand_numbers=np.random.rand((100000))
            rand_numbers2=np.random.rand((100000))
            # sample the sizes of the impactors
            diams, mass = impactor_model.sample_sizes(rand_numbers, rho_pr)
            # sample the impact velocities - not sure if this should be a different random number
            vels = impactor_velocity.inv_cdf1(rand_numbers2)
        
        
        # solve ODEs
        #y0 = ystore[i,0:2]
        #sol = odeint(test, y0, [t[i], t[i+1]], args=(omega,) )
        
        # 1. collapse? 0.5 bar and assumed to be 6 mbar
        
        # 2. calculate H2O - two temperatures
        
        # 3. impacts of asteroids and comets. during this time-step
        num_impactor1=crater_chronology.N1_between(Amars,-t[i+1]/1e9,-t[i]/1e9)
        num_impactor20=crater_chronology.N20_between(Amars,-t[i+1]/1e9,-t[i]/1e9,model=1)
        

        # note that results at each time get scaled by:
        # result * total_impactor_mass / np.sum(mass) * num_impactor1 / num_impactor1_tot
        scaling = total_impactor_mass / np.sum(mass) * num_impactor1 / num_impactor1_tot
        
        
        
        # 4. loss of the atmosphere following empirical formulae
        # sum(Matm) / sum(Natm) is the molecular weight
        rho_0 = Patm*1.e5/tmars/(Rgas/(np.sum(Matm)/np.sum(Natm)))
        xi=np.maximum(impactor_atmos_loss_gain_eq4_5. \
            dimensionless_xi(diams,rho_t,rho_pr, vels, uesc, Hscale, rho_0)  , 1.e-10) 
        deltaM=impactor_atmos_loss_gain_eq4_5.change_atmos(xi,vels,uesc,mass)
        
        
        # left in the projectiles is (1-Xpr)*mass, but then only a certain % of this 
        #Xpr=np.maximum(impactor_atmos_loss_gain_eq4_5. \
        #    normalised_projectile_mass(rho_t,rho_pr,vels,uesc,xi),0)
        #deltaM=deltaM-Xpr*mass*0.01
        
        frac = Matm / np.sum(Matm) # initial fraction of each component - stays the same
        Matm_final = np.sum(Matm) - np.sum(deltaM)*scaling
        Matm = frac *Matm_final
        Natm = Matm /MolW_atm
        
        # 5. Sputtering and Photochemical escape
        (euv_flux, f_co2_flux) = sputtering_co2.sputtering_co2_rate(-t[i]/1e9)
        Natm[0] = Natm[0] - f_co2_flux * tstep *86400*365/6.02e23
        Natm[0] = np.maximum(Natm[0],0.)
        
        # 6. Volcanic degassing - digitise rates and incorporate total
        Matm = Natm * MolW_atm
        h2o=volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'h2o')
        co2=volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'co2')
        n2 =volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'n2')
        Matm = Matm + np.array([co2,n2,h2o])*MolW_atm*tstep*86400*365/6.02e23
        
        # 7. Isotope fractionation. 
        
        
        
        # readjust Matm
        Patm = np.maximum(np.sum(Matm)*grav_mars/Amars/1.e5, 6e-3)
        # if Patm < 0.5 bar - the CO2 condenses at the poles - only the CO2 so will have 
        # to wait to do properly
        if(Patm < Pcollapse):
            Patm = 6.e-3
            Matm[0] = Patm*1e5*Amars / grav_mars # mass of atmosphere
            Matm[1:] = 0. # if it's collapsed, just assume it is CO2

        
        
        Natm = Matm / MolW_atm
        
        # update scale height
        Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
        
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


