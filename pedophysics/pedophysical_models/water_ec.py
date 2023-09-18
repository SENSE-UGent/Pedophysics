def SenGoode(T, C_f):
    """
    Calculate soil water real electrical conductivity using the Sen and Goode model.

    Based on the model proposed by Sen and Goode (1992 a, b) [1][2], this function estimates 
    the oil water real electrical conductivity given the soil temperature and 
    salinity of the bulk pore fluid.

    Parameters
    ----------
    T : array_like
        Soil bulk temperature [K]
    C_f : array_like
        Soil (NaCl) salinity of the bulk pore fluid [mol/L].

    Returns
    -------
    water_ec : array_like
        Soil water real electrical conductivity [S/m]].

    References
    ----------
    .. [1] Sen P and Goode P (1992a) 
    Influence of temperature on electrical conductivity on shaly sands. 
    Geophysics 57: 89-96.
    .. [2] Sen P and Goode P (1992b) 
    Errata, to: Influence of temperature on electrical conductivity of shaly sands. 
    Geophysics 57: 1658.

    Notes
    -----
    The function uses specific coefficients based on the Sen and Goode model to 
    compute the electrical conductivity of soil water, taking into account the 
    temperature and salinity effects.

    Example
    -------
    >>> SenGoode(298.15, 0.01)
    0.117822
    
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