import numpy as np

IDP_table3 = dict()
IDP_table3['Ne'] = 4.89e12
IDP_table3['Ar'] = 1.11e12
IDP_table3['Kr'] = 1.64e10
IDP_table3['Xe'] = 7.33e9

IDP_ar = dict()
IDP_ar['Ne'] = 20
IDP_ar['Ar'] = 36
IDP_ar['Kr'] = 84
IDP_ar['Xe'] = 130



def get_mole(elements,Ne=1):
    rates = np.zeros(len(elements))
    
    for i in range(len(elements)):
        rates[i] = IDP_table3[elements[i]] / IDP_ar[elements[i]]
        if(elements[i] == 'Ne'):
            rates[i] *= Ne
            
    return rates