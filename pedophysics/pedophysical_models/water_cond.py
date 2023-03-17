def SenGoode(T, c):
    """
        Sen and Goode (1992, a, b)
        
        Parameters
        ----------
        temperature: float
            Soil temperature [K]
        
        salinity: float
            Soil (NaCl) salinity of the bulk pore fluid [mol/L]
 
        Returns
        -------
        water_cond: float
            Soil water real electrical conductivity [S/m]   
    """   
    T_celsius = T-273.15 
    d1 = 5.6 
    d2 = 0.27 
    d3 = -1.51e-4 
    d4 = 2.36 
    d5 = 0.099 
    d6 = 0.214 
    water_cond = (d1+d2*T_celsius+d3*T_celsius**2)*c - ((d4+d5*T_celsius)/(1+d6*c**0.5))*c**1.5

    return water_cond