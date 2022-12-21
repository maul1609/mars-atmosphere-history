import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import odeint

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
    tfinal=0.
    tstep=1e6
    omega=2.*np.pi*1.e-9
    output_interval = 1e9
    last_output = tinit-output_interval
    """
        ----------------------------------------------------------------------------------
    """

    """
        auxiliary calcs ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    t = np.mgrid[tinit:tfinal+tstep:tstep]
    nsteps = len(t)
    n_ode=3
    """
        ----------------------------------------------------------------------------------
    """

    """
        set-up arrays ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    ystore=np.zeros((nsteps,n_ode))
    ystore[0,0]=1.
    ystore[0,1]=0.
    """
        ----------------------------------------------------------------------------------
    """


    """
        main loop ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    """
    for i in range(nsteps-1):
        if ((t[i]-last_output)>=output_interval):
            print('time: ' + str(t[i]/1e9))
            last_output = t[i]
                
        y0 = ystore[i,0:2]
        
        
        # solve ODEs
        sol = odeint(test, y0, [t[i], t[i+1]], args=(omega,) )
        
        # 1. collapse? 0.5 bar and assumed to be 6 mbar
        
        # 2. calculate H2O - two temperatures
        
        # 3. impacts of asteroids and comets. during this time-step
        
        # 4. loss of the atmosphere following empirical formulae
        
        # 5. Sputtering and Photochemical escape
        
        # 6. Volcanic degassing - see rates
        
        # 7. Isotope fractionation. 
        
        
        
        
        ystore[i+1,0]=sol[-1,0]
        ystore[i+1,1]=sol[-1,1]
        
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


