# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 12:02:03 2022

@author: purpl
"""
import numpy as np

def F_i_sp(F_CO2,Y_i,Y_CO2,N_i,N_CO2,R_diff,alpha):
    """
    See equation 6 from Kurokawa et al (2018)
    
    i  - the gas species 
    sp - sputtering rate
    F_i_sp - sputtering rate of that gas species
    Y_i - the yield 
    N_i - the number density at homopause
    R_diff,_i_CO2 - the fractionatin factor by diffusive separation between 
    homopause and exobase
    alpha - dimensionlass factor, defined by 
    alpha = sumof_i(N_i/N_CO2,R_diff,_i/_CO2)
    
    alpha denotes the summation of the mixing ratios at the exobase. 
    The fractionation by diffusive separation is given by
    equation 7

    """
    #the number of commas in this equation are confusing me as 
    #I'm unsure what they're denoting.
    
    F_i_sp1=(F_CO2*(Y_i/Y_CO2)*(N_i/N_CO2)*(R_diff)*(1/alpha))
    
    return F_i_sp1

def alpha(N_i,N_CO2,R_diff):
    
    alpha1=np.sum((N_i/N_CO2)*(R_diff))
    
    return alpha1

def R_diffij(exp,deltam_i_j,g,deltaz,k_B,T):
    
    #I'm unsure here as R_diff is written using the R_diff,i to denote species
    # and R_diff,i/j, how do I put that in def due to the use of operators? 
    
    """
    j           - CO2 for R_diff,i/CO2
    deltam_i_j  - is the mass difference between two species
    """
    
    R_diff=np.exp(-deltam_i_j*g*deltaz/(k_B*T))
    
    return R_diff



#table 1, Fco2, for equation 6, in Luhmann et al 1992. 