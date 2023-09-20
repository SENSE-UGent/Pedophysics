def MalmbergMaryott(T):
    """
    Calculate soil water phase real dielectric permittivity using the Malmberg & Maryott model.

    This function utilizes the model proposed by Malmberg and Maryott (1956) [1] to estimate 
    the soil water phase real dielectric permittivity based on a given soil temperature. Reported RMSE = 0.0046

    Parameters
    ----------
    T : array_like
        Soil bulk temperature [K].

    Returns
    -------
    water_perm : array_like
        Soil water phase real dielectric permittivity [-]

    References
    ----------
    .. [1] Malmberg C and Maryott A (1956) 
    Dielectric constant of water from 0 °C to 100 °C. 
    Journal of Research of the National Bureau of Standards 56(1): 1-8, Paper 2641.

    Example
    -------
    >>> MalmbergMaryott(298.15)
    78.303
    """
    T_c = T - 273.15 # Kelvin to Celsius
    water_perm = 87.740 - 0.40008*T_c + 9.398e-4*T_c**2 - 1.410e-6*T_c**3 

    return water_perm


def Olhoeft(T, C_f):
    """
    Calculate soil water phase real dielectric permittivity using the Olhoeft (1986) [1] model.

    Parameters
    ----------
    T : array_like
        Soil bulk temperature [K].
    C_f : array_like
        Soil (NaCl) salinity of the bulk pore fluid [mol/L].

    Returns
    -------
    water_perm : array_like
        Soil water phase real dielectric permittivity [-]

    References
    ----------
    .. [1] Revil, A., Schwaeger, H., Cathles, L. M., and Manhardt, P. D. (1999), 
    Streaming potential in porous media: 2. Theory and application to geothermal systems, 
    J. Geophys. Res., 104(B9), 20033-20048, doi:10.1029/1999JB900090.

    Example
    -------
    >>> Olhoeft(298.15, 0.1)
    76.945
    """
    a0 = 295.68
    a1 = -1.2283 
    a2 = 2.094e-3
    a3 = -1.41e-6
    c1 = -13 
    c2 = 1.065
    c3 = -0.03006 
    water_perm = a0 + a1*T + a2*T**2 + a3*T**3 + c1*C_f + c2*C_f**2 + c3*C_f**3

    return water_perm
    