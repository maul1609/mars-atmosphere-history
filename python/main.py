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
import MarsSputtering_Calc6_7
import isotopic_data
import photo_chemical
import IDPs

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
    tmars = 270. # temperature of martian surface - may change???
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
    
    # IDP rates, table 2
    IDP_rate_moles = IDPs.get_mole(['Ne','Ar','Kr','Xe'],Ne=1)
    
    isotope_comps_sim = np.array(isotopic_data.isotope_comps['Volcanic degassing'])
    ind,=np.where(['Xe' in iso for iso in isotopic_data.isotopes])
    isotope_comps_sim[ind]=np.array(isotopic_data.isotope_comps['Mars'])[ind]
    """
        ----------------------------------------------------------------------------------
    """




    """
        auxiliary calcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    Amars=4.*np.pi*Rmars**2 # surface area of mars
    Matm = Patm*1e5*Amars / grav_mars # mass of atmosphere
    # call function to set isotopes and adjust atmosphere to same pressure
    (mole1,mole_elements)=isotopic_data.set_amount_of_element(Matm/MolW_atm[0])
    Natm = np.concatenate([mole1,[0]]) # number of moles of co2, n2, h2o
    Matm = Natm *MolW_atm  # mass of co2, n2, h2o
#     Natm = Matm / MolW_atm     # number of moles of co2, n2, h2o
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
        dM = np.sum(deltaM)*scaling
        const1 = dM / np.sum(MolW_atm*Natm)
        Natm = Natm - const1*Natm
        Matm = Natm*MolW_atm
        # remove elements in the impact too
        mole_elements = mole_elements*(1.-const1)
        
        
        # 5. Sputtering and Photochemical escape
        # sputter each element / isotope: c (Co2), N (N2), Ne, Ar, Kr, Xe       
        (euv_flux, f_co2_flux) = sputtering_co2.sputtering_co2_rate(-t[i]/1e9)
        F_i_sp1=MarsSputtering_Calc6_7.F_i_sp(f_co2_flux, \
            mole_elements,Natm[0],['C','N','Ne','Ar','Kr','Xe'])
        # CO2 and N2
        Natm[0] = Natm[0] - F_i_sp1[0] * tstep *86400*365/6.02e23
        Natm[1] = Natm[1] - 0.5*F_i_sp1[1] * tstep *86400*365/6.02e23
        Natm = np.maximum(Natm,0.)
        # now elements
        mole_elements[0] = mole_elements[0] - F_i_sp1[0] * tstep *86400*365/6.02e23
        mole_elements[1] = mole_elements[1] - F_i_sp1[1] * tstep *86400*365/6.02e23
        mole_elements[2] = mole_elements[2] - F_i_sp1[2] * tstep *86400*365/6.02e23
        mole_elements[3] = mole_elements[3] - F_i_sp1[3] * tstep *86400*365/6.02e23
        mole_elements[4] = mole_elements[4] - F_i_sp1[4] * tstep *86400*365/6.02e23
        mole_elements[5] = mole_elements[5] - F_i_sp1[5] * tstep *86400*365/6.02e23
        
        
        # 5b photochemical escape
        # carbon - rate of escape
        Fcph = photo_chemical.carbon_escape((t[i]-tinit)/1e9,\
            (t[i+1]-tinit)/1e9,(-tinit)/1e9) / 6.02e23
        # nitrogen - rate of escape
        (fn2_1,fn2_2,fn2_3)=photo_chemical.nitrogen_escape(np.maximum(Natm[0],1e-3),Natm[1],euv_flux,Amars)
        # total loss rate
        fn2 = fn2_1+fn2_2+fn2_3
        Natm[0] = Natm[0] - Fcph*tstep*86400*365/6.02e23
        Natm[1] = Natm[1] - fn2*tstep*86400*365/6.02e23
        mole_elements[0] -= Fcph*tstep*86400*365/6.02e23
        mole_elements[1] -= fn2*tstep*86400*365/6.02e23 * 2
        Matm = Natm * MolW_atm
        
        
        
        # 5c IDPs - just table 3, Kurokawa et al.
        mole_elements[2:] += IDP_rate_moles*tstep / 1e9
        
        
        
        
        # 6. Volcanic degassing - digitise rates and incorporate total
        Matm = Natm * MolW_atm
        h2o=volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'h2o')
        co2=volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'co2')
        n2 =volcano_outgassing.get_outgas_rate(t[i+1]/1e9,'n2')
        Matm = Matm + np.array([co2,n2,h2o])*MolW_atm*tstep*86400*365/6.02e23
        Matm[2] = np.minimum(Matm[2]*grav_mars/Amars,ph2o)*Amars/grav_mars
        Natm = Matm / MolW_atm
        # 7. Isotope fractionation. 
        
        
        
        # readjust Matm
        Patm = np.maximum(np.sum(Matm)*grav_mars/Amars/1.e5, 6e-3)
        # if Patm < 0.5 bar - the CO2 condenses at the poles - only the CO2 so will have 
        # to wait to do properly. e.g. need to know how lack of absorption by CO2 will affect 
        # temperature, and then link this to the clausius clapeyron equation
        # same for water vapour really, perhaps having a reservoir
        if(Patm < Pcollapse):
            Patm = 6.e-3
            Matm[0] = Patm*1e5*Amars / grav_mars # mass of atmosphere
            Patm = np.sum(Matm)*grav_mars/(Amars)/1e5
            #Matm[1:] = 0. # if it's collapsed, just assume it is CO2
            mole_elements[0] = Matm[0] / MolW_atm[0]

        
        
        Natm = Matm / MolW_atm
        mole_elements[0] += co2*tstep*86400*365/6.02e23
        mole_elements[1] += 2*n2*tstep*86400*365/6.02e23
        
        # update scale height
        Hscale=Rgas*tmars/(np.sum(Matm)/np.sum(Natm)*grav_mars) # scale height - m - may change???
        
        #ystore[i+1,0]=sol[-1,0]
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


