def SenGoode(T, C_f):
    """
        Sen and Goode (1992 a, b)
        
        Parameters
        ----------
        temperature: float
            Soil temperature [K]
        
        C_f: float
            Soil salinity (NaCl) of the bulk pore fluid [mol/L]
 
        Returns
        -------
        water_ec: float
            Soil water real electrical conductivity [S/m]   
    """   
    T_celsius = T-273.15 
    d1 = 5.6 
    d2 = 0.27 
    d3 = -1.51e-4 
    d4 = 2.36 
    d5 = 0.099 
    d6 = 0.214 
    water_ec = (d1+d2*T_celsius+d3*T_celsius**2)*C_f - ((d4+d5*T_celsius)/(1+d6*C_f**0.5))*C_f**1.5

    return water_ec


#def Hilhorst(bulk_ec, bulk_perm, wp, bulk_perm_0):
    """
        Hilhorst 2000
        
        Parameters
        ----------
        bulk_ec: float
            Soil bulk real electrical conductivity [S/m]

        bulk_perm: float
            Soil bulk real relative dielectric permittivity   
        
        wp: float
            water permittivity phase [-]

        bulk_perm_0: float
            Soil bulk real relative dielectric permittivityat zero bulk real electrical conductivity [-]   
 
        Returns
        -------
        water_ec: float
            Soil water real electrical conductivity [S/m]   
    """   
#    water_ec = 54444 # not sure of I should write this
#    return water_ec