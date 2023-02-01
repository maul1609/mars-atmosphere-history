import numpy as np

# solar abundances of these isotopes, table 6, Lodder, 20xx
solar_abundances=np.array([7.001e6, 1.943e6, 1.996e6, 86710, 31.38, 0.236])
carbon_12_pc = 98.8922
nitrogen_14_pc = 99.6337
solar_percents = np.array([98.8922, 99.6337, 92.9431, 85.5946, 56.90, 4.38])


# figure 2, Kurokawa et al
abundances1 = dict()
volx = np.array([12, 14, 20, 36, 84, 130])
voly = np.array([2.257863844,2.757493323,0.000446188,0.949493424,55.30511164,16.66728846])
abundances1['Comets'] = np.array([volx,voly])

volx = np.array([12,14, 20, 36, 84, 130])
voly = np.array([np.nan, np.nan, 0.005999776,0.002360426,0.057817309,1.729585844])
abundances1['IDPs'] = np.array([volx,voly])

volx = np.array([12, 14, 20, 36, 84, 130])
voly = np.array([0.727338094,0.120333307,1.36E-08,1.03771E-06,4.63015E-05,0.000665505])
abundances1['Asteroids'] = np.array([volx,voly])

voly = np.array([3.0137E-06,4.66E-07,2.97E-11,3.37E-09,3.13E-07,2.56E-07])
abundances1['Mars'] = np.array([volx,voly])

voly = np.array([0.000175539, 8.1881E-06, 4.14E-11,3.88E-11,7.49E-09,1.56E-08])
abundances1['Volcanic degassing'] = np.array([volx,voly])

atom_wgts=volx/1000.



del volx,voly

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
            
isotopes[1] = '22Ne/20Ne'
isotopes[2] = '38Ar/36Ar'
isotopes[3] = '78Kr/80Kr'
isotopes[4] = '82Kr/80Kr'
isotopes[5] = '83Kr/80Kr'
isotopes[6] = '84Kr/80Kr'
isotopes[7] = '86Kr/80Kr'

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


isotope_comps['Mars'][1]=1./isotope_comps['Mars'][1]
isotope_comps['Mars'][2]=1./isotope_comps['Mars'][2]
rat1 = isotope_comps['Mars'][3] # 78/84
rat2 = isotope_comps['Mars'][4] # 80/84
rat3 = isotope_comps['Mars'][5] # 82 /84
rat4 = isotope_comps['Mars'][6] # 83 /84
rat4 = isotope_comps['Mars'][7] # 86 /84
isotope_comps['Mars'][3]=rat1 / rat2 # 78 /80
isotope_comps['Mars'][4]=rat3 / rat2 # 82 /80
isotope_comps['Mars'][5]=rat4 / rat2 # 83 /80
isotope_comps['Mars'][6]=1./rat2 # 84/80
isotope_comps['Mars'][7]=rat4 / rat2 # 86/80

isotope_comps['Volcanic degassing'][1]=1./isotope_comps['Volcanic degassing'][1]
isotope_comps['Volcanic degassing'][2]=1./isotope_comps['Volcanic degassing'][2]
rat1 = isotope_comps['Volcanic degassing'][3] # 78/84
rat2 = isotope_comps['Volcanic degassing'][4] # 80/84
rat3 = isotope_comps['Volcanic degassing'][5] # 82 /84
rat4 = isotope_comps['Volcanic degassing'][6] # 83 /84
rat4 = isotope_comps['Volcanic degassing'][7] # 86 /84
isotope_comps['Volcanic degassing'][3]=rat1 / rat2 # 78 /80
isotope_comps['Volcanic degassing'][4]=rat3 / rat2 # 82 /80
isotope_comps['Volcanic degassing'][5]=rat4 / rat2 # 83 /80
isotope_comps['Volcanic degassing'][6]=1./rat2 # 84/80
isotope_comps['Volcanic degassing'][7]=rat4 / rat2 # 86/80

isotope_comps['Asteroids'][1]=1./isotope_comps['Asteroids'][1]
isotope_comps['Asteroids'][2]=1./isotope_comps['Asteroids'][2]
rat1 = isotope_comps['Asteroids'][3] # 78/84
rat2 = isotope_comps['Asteroids'][4] # 80/84
rat3 = isotope_comps['Asteroids'][5] # 82 /84
rat4 = isotope_comps['Asteroids'][6] # 83 /84
rat4 = isotope_comps['Asteroids'][7] # 86 /84
isotope_comps['Asteroids'][3]=rat1 / rat2 # 78 /80
isotope_comps['Asteroids'][4]=rat3 / rat2 # 82 /80
isotope_comps['Asteroids'][5]=rat4 / rat2 # 83 /80
isotope_comps['Asteroids'][6]=1./rat2 # 84/80
isotope_comps['Asteroids'][7]=rat4 / rat2 # 86/80

isotope_comps['Comets'][1]=1./isotope_comps['Comets'][1]
isotope_comps['Comets'][2]=1./isotope_comps['Comets'][2]
rat1 = isotope_comps['Comets'][3] # 78/84
rat2 = isotope_comps['Comets'][4] # 80/84
rat3 = isotope_comps['Comets'][5] # 82 /84
rat4 = isotope_comps['Comets'][6] # 83 /84
rat4 = isotope_comps['Comets'][7] # 86 /84
isotope_comps['Comets'][3]=rat1 / rat2 # 78 /80
isotope_comps['Comets'][4]=rat3 / rat2 # 82 /80
isotope_comps['Comets'][5]=rat4 / rat2 # 83 /80
isotope_comps['Comets'][6]=1./rat2 # 84/80
isotope_comps['Comets'][7]=rat4 / rat2 # 86/80

isotope_comps['IDPs'][1]=1./isotope_comps['IDPs'][1]
isotope_comps['IDPs'][2]=1./isotope_comps['IDPs'][2]
rat1 = isotope_comps['IDPs'][3] # 78/84
rat2 = isotope_comps['IDPs'][4] # 80/84
rat3 = isotope_comps['IDPs'][5] # 82 /84
rat4 = isotope_comps['IDPs'][6] # 83 /84
rat4 = isotope_comps['IDPs'][7] # 86 /84
isotope_comps['IDPs'][3]=rat1 / rat2 # 78 /80
isotope_comps['IDPs'][4]=rat3 / rat2 # 82 /80
isotope_comps['IDPs'][5]=rat4 / rat2 # 83 /80
isotope_comps['IDPs'][6]=1./rat2 # 84/80
isotope_comps['IDPs'][7]=rat4 / rat2 # 86/80



def set_amount_of_element(carbon_moles):
    """
        -from moles of total carbon, work out moles of carbon-12
        -from moles of carbon-12 you can calculate moles of others
        -and convert to mass using atomic weights
    """
    
    n_moles = np.zeros(len(atom_wgts))
    # carbon-12
    n_moles[0] = carbon_moles * carbon_12_pc/100.
    for i in range(1,len(atom_wgts)-1):
        # N, Ne, Ar, Kr
        n_moles[i] = carbon_moles * carbon_12_pc/100. * \
            abundances1['Volcanic degassing'][1][i] / abundances1['Volcanic degassing'][1][0] * \
            solar_abundances[i] / solar_abundances[0]
    
    # Xe
    n_moles[len(atom_wgts)-1] = carbon_moles * carbon_12_pc/100. * \
        abundances1['Mars'][1][len(atom_wgts)-1] / abundances1['Mars'][1][0] * \
        solar_abundances[len(atom_wgts)-1] / solar_abundances[0]


    # molecular nitrogen moles
    n2_moles = (n_moles[1]/2.) / (nitrogen_14_pc/100.)
    
    # scale all moles so that they equal the initial number
    initial_number = carbon_moles
    
    denominator = (carbon_moles+n2_moles+np.sum(n_moles[2:]))
    # adjusted carbon
    carbon_moles = carbon_moles * initial_number / denominator
    # adjusted nitrogen
    n2_moles = n2_moles * initial_number / denominator
    # adjusted everything else
    n_moles = n_moles * initial_number / denominator
    
    
    
    
    # now return the total mass of carbon dioxide and molecular nitrogen
    # and also return the total number of moles of the chosen isotopes 
    
    return ([carbon_moles,n2_moles],n_moles)
    

def continuous_sources(Imn, Imn0, N, dN):
    Imn_new=[None]*len(Imn)
    # equation 10
    for m in range(len(Imn)):
        Imn_new[m] = [None] * len(Imn[m])
        for n in range(len(Imn[m])):
            Imn_new[m][n]=Imn[m][n]+(Imn0[m][n]-Imn[m][n])*np.sum(Imn[m]) / \
                np.sum(Imn0[m]) * dN[m] / N[m]

    return Imn_new
    


    
def impact_replenishment(Imn, Imn0_asteroids, Imn0_comets, N, dN, dN2):
    Imn_new=[None]*len(Imn)
    # asteroids eqaution 11
    for m in range(len(Imn)):
        Imn_new[m] = [None] * len(Imn[m])
        for n in range(len(Imn[m])):
            Imn_new[m][n]=(Imn[m][n]*np.sum(Imn0_asteroids[m])*N[m] + \
                Imn0_asteroids[m][n]*np.sum(Imn[m])*dN[m]) / \
                (np.sum(Imn0_asteroids[m])*N[m] + np.sum(Imn[m])*dN[m])

    Imn_new2=[None]*len(Imn)
    # comets
    N += dN
    for m in range(len(Imn)):
        Imn_new2[m] = [None] * len(Imn[m])
        for n in range(len(Imn[m])):
            Imn_new2[m][n]=(Imn_new[m][n]*np.sum(Imn0_comets[m])*N[m] + \
                Imn0_comets[m][n]*np.sum(Imn_new[m])*dN2[m]) / \
                (np.sum(Imn0_comets[m])*N[m] + np.sum(Imn_new[m])*dN2[m])

    return Imn_new2
    
if __name__ == "__main__":
    """
        plot
    """
    
    
    import matplotlib.pyplot as plt
    plt.ion()
    plt.plot(abundances1['Comets'][0],abundances1['Comets'][1])
    plt.plot(abundances1['IDPs'][0],abundances1['IDPs'][1])
    plt.plot(abundances1['Asteroids'][0],abundances1['Asteroids'][1])
    plt.plot(abundances1['Mars'][0],abundances1['Mars'][1])
    plt.plot(abundances1['Volcanic degassing'][0],abundances1['Volcanic degassing'][1])

    plt.legend(abundances1.keys())
    plt.yscale('log')   
    plt.xticks([12,14,20,36,84,130])    
    plt.gca().set_xticklabels(['C','N','Ne','Ar','Kr','Xe'])  