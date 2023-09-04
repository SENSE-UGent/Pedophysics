def MalmbergMaryott(T):
    """
        Malmberg & Maryott, 1956. RMSE = 0.0046
        
        Parameters
        ----------
        temperature: float
            Soil temperature [K]
 
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]
    """    
    T_c = T - 273.15 # Kelvin to Celsius
    water_perm = 87.740 - 0.40008*T_c + 9.398e-4*T_c**2 - 1.410e-6*T_c**3

    return water_perm


def Olhoeft(T, c):
    """
        Olhoeft, 1986
        
        Parameters
        ----------
        temperature: float
            Soil temperature [K]
        
        salinity: float
            Soil (NaCl) salinity of the bulk pore fluid [mol/L]
 
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]
    """    
    a0 = 295.68
    a1 = -1.2283 
    a2 = 2.094e-3
    a3 = -1.41e-6
    c1 = -13 
    c2 = 1.065
    c3 = -0.03006 
    water_perm = a0 + a1*T + a2*T**2 + a3*T**3 + c1*c + c2*c**2 + c3*c**3

    return water_perm


def Stogryn(T, c, f):
    """
        Stogryn, 1971
        
        Parameters
        ----------
        temperature: float
            Soil temperature [K]
        
        salinity: float
            Soil (NaCl) salinity of the bulk pore fluid [mol/L]
 
        frequency_perm: float
            frequency of permittivity readings [Hz]
            
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]
    """
    T_c = T - 273.15 # Kelvin to Celsius
    print(T_c)
    molarity = 58.44 # Molar mass of NaCl
    water_perm_inf = 4.5

    c_ppt = c*molarity

    print(c_ppt)

    N = c_ppt*(1.707e-2 + 1.205e-5*c_ppt + 4.058e-9*c_ppt**2)
    a_N = 1 - 0.2551*N + 5.151e-2*N**2 - 6.889e-3*N**3
    b_N_T_c = 0.1463e-2*N*T_c + 1 - 0.04896*N - 0.02967*N**2 + 5.644e-3*N**3
    e_t_0 = 87.74 + 4.0008*T_c + 9.398e-4*T_c**2 + 1.41e-6*T_c**3
    two_pi_tau_T_c_0 = 1.1109e-10 - 3.824e-12*T_c + 6.938e-14*T_c**2 - 5.096e-16*T_c**3
    water_perm = water_perm_inf + (e_t_0*a_N - water_perm_inf)/(1 + (two_pi_tau_T_c_0*b_N_T_c*f)**2)

    return water_perm
    