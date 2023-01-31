import numpy as np
import scipy.interpolate as scint 

    
def sputtering_co2_rate(timeGyr):
    """
        See Luhmann et al. 1992, table 1
        calculates the sputtering rate of CO2 in Martian atmosphere
    """
    euv_grid = np.array([6., 3., 1.])
    time_grid = np.array([1., 2., 4.5])
    f_co2_grid = np.array([3.e26, 6.e25, 3.e25])
    
    # generate my interpolation function from the data points
    euv_func=scint.interp1d(time_grid,euv_grid,fill_value='extrapolate',kind='linear')
    f_co2_func=scint.interp1d(time_grid,f_co2_grid,fill_value='extrapolate',kind='linear')
    
    
    
    return (euv_func(timeGyr), f_co2_func(timeGyr) )
    
    
if __name__ == "__main__":
    # test the function by evaluating and plotting
    import matplotlib.pyplot as plt 
    
    t_test_gyr = np.linspace(0.,5.,100)
    
    (euv_flux, f_co2_flux) = sputtering_co2_rate(t_test_gyr)
    
    plt.ion()
    plt.figure()
    plt.subplot(211)
    plt.plot(t_test_gyr,f_co2_flux)    
    plt.xlabel('time (Gyr)')
    plt.ylabel('F$_{CO2}$')
    plt.subplot(212)
    plt.plot(t_test_gyr,euv_flux)
    plt.xlabel('time (Gyr)')
    plt.ylabel('Relative EUV flux')
    plt.show()