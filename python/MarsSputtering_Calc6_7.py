# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:02:03 2022

@author: purpl
"""
import numpy as np

k_B=1.381e-23 # Boltzmann constant
R_gas = 8.314 # ideal gas constant
Navo  = 6.02e23

Y_CO2=0.67 # see page 3 of Luhmann et al.

amu=1.66e-27
g_mars=3.721
deltaz=40.*1000. # 40000 m
T=100.
amu_co2=44

yields = dict()
yields['C'] = 0.67 # strictly 0.7
yields['N'] = 2.4
yields['Ne'] = 3.
yields['Ar'] = 1.4
yields['Kr'] = 1.143
yields['Xe'] = 0.738

amus = dict()
amus['C'] = 12
amus['N'] = 14
amus['Ne'] = 20
amus['Ar'] = 36
amus['Kr'] = 84
amus['Xe'] = 130

def calculate_rdiffs(element_mass,isotope_mass):
    """
       easy way of calculating Rdiffs factors 
    """
    Rdiffs = [None]*len(element_mass)
    for i in range(len(element_mass)):
        Rdiffs[i] = [None]*len(isotope_mass[i])
        for j in range(len(isotope_mass[i])):
            delta_m=(element_mass[i]-isotope_mass[i][j])*amu
            Rdiffs[i][j] = R_diffij(delta_m,g_mars,deltaz,T)
            
    return Rdiffs
    
def F_i_sp(F_CO2,N_i,N_CO2,elements):
    """
    See equation 6 from Kurokawa et al (2018)
    
    i  - the gas species 
    sp - sputtering rate
    F_i_sp - sputtering rate of that gas species
    Y_i - the yield 
    N_i - the number density at homopause
    R_diff,_i_CO2 - the fractionatin factor by diffusive separation between 
    homopause and exobase
    alpha1 - dimensionlass factor, defined by 
    alpha1 = sumof_i(N_i/N_CO2,R_diff,_i/_CO2)
    
    alpha1 denotes the summation of the mixing ratios at the exobase. 
    The fractionation by diffusive separation is given by
    equation 7

    """
    #the number of commas in this equation are confusing me as 
    #I'm unsure what they're denoting.
    Y_i = np.zeros(len(elements))
    delta_m = np.zeros(len(elements))
    R_diff = np.zeros(len(elements))
    for i in range(len(elements)):
        Y_i[i] = yields[elements[i]]
        delta_m[i]=(amus[elements[i]]-amu_co2)*amu
        R_diff[i] = R_diffij(delta_m[i],g_mars,deltaz,T)
    
    
    alpha1 = alpha(N_i, N_CO2, R_diff)
    
    F_i_sp1=(F_CO2*(Y_i/Y_CO2)*(N_i/N_CO2)*(R_diff)*(1./alpha1))
    
    return F_i_sp1

def alpha(N_i,N_CO2,R_diff):
    
    alpha1=np.sum((N_i/N_CO2)*(R_diff))
    
    return alpha1

def R_diffij(deltam_i_j,g,deltaz,T):
    
    #I'm unsure here as R_diff is written using the R_diff,i to denote species
    # and R_diff,i/j, how do I put that in def due to the use of operators? 
    
    """
    j           - CO2 for R_diff,i/CO2
    deltam_i_j  - is the mass difference between two species
    """
    
    R_diff=np.exp(-deltam_i_j*g*deltaz/(k_B*T))
    
    return R_diff



#table 1, Fco2, for equation 6, in Luhmann et al 1992. 

if __name__ == "__main__":
    
    # Ideal gas law
    # P = (N/V) R T
    # (N/V) = P / (R T)
    Press = 0.5e5
    Temp = 273.15
    Ntot = Press / (R_gas * Temp) * Navo
    
    print(Ntot)
    
    
    # test the fractionation by diffusive separation 
    # tested against table 3 in Jakosky et al. 1994
    # watch out for whether it is a molecule or atom
    amu_n2=28
    
    amu_n_15=15
    amu_n_14=14

    # R(28/44)
    delta_m=(amu_n2-amu_co2)*amu
    R_n2_co2=R_diffij(delta_m,g_mars,deltaz,T)

    # R(15/14)
    delta_m=(amu_n_15-amu_n_14)*amu
    R_n_14_n_15=R_diffij(delta_m,g_mars,deltaz,T)
    
    # R(13/12)
    amu_c_13=13
    amu_c_12=12
    delta_m=(amu_c_13-amu_c_12)*amu
    R_c_13_c_12=R_diffij(delta_m,g_mars,deltaz,T)

    
    print(R_n2_co2, R_n_14_n_15, R_c_13_c_12)
    
    
    # alpha calculation
    amu_species = [16, 20, 40, 44, 28]
    species = ['O', 'Ne', 'Ar', 'CO2', 'N2']
    fracall = np.array([0.01, 0.01, 0.01, 0.96, 0.01])
    # calculate RdiffiCO2
    RdiffiCO2 = np.zeros(len(fracall))
    for i in range(len(fracall)):
        delta_m = amu*(amu_species[i]-amu_species[3])
        RdiffiCO2[i]=R_diffij(delta_m,g_mars,deltaz,T)
        
    print(RdiffiCO2)
    
    # alpha can now be calculated
    Ni = fracall * Ntot
    alpha1 = alpha(Ni, Ni[3], RdiffiCO2)
    print(alpha1)
    
    