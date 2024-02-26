import numpy as np

def WunderlichP(water, perm_init, wat_init, wp, Lw): 
    """
    Calculate the soil bulk real relative dielectric permittivity using the Wunderlich model.

    This is a effective medium model that uses a differential 
    approach to compute the dielectric permittivity based on the initial 
    conditions and the water content [1]. Reported RMSE = 0.009

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    perm_init : float
        Initial soil bulk real relative dielectric permittivity [-].
    wat_init : float
        Initial soil volumetric water content [m**3/m**3].
    wp : array_like
        Soil water phase real dielectric permittivity [-].
    Lw : float
        Soil scalar depolarization factor of water aggregates (effective medium theory) [-]

    Returns
    -------
    array_like
        The estimated soil bulk real relative dielectric permittivity [-]

    Notes
    -----
    The Wunderlich model is a differential model and the approach used here 
    employs a simple while loop with a step of 0.01 until the differential 
    fraction reaches 1. The method takes into account the initial water content, 
    initial relative dielectric permittivity, and weighting factors to determine the 
    bulk real relative dielectric permittivity.

    References
    ----------
    .. [1] Tina Wunderlich, Hauke Petersen, Said Attia al Hagrey, Wolfgang Rabbel; 
    Pedophysical Models for Resistivity and Permittivity of Partially Water-Saturated Soils. 
    Vadose Zone Journal 2013;; 12 (4): vzj2013.01.0023. doi: https://doi.org/10.2136/vzj2013.01.0023

    Example
    -------
    >>> WunderlichP(0.3, 7, 0.05, 80, 0.01)
    24.591

    """
    diff = water - wat_init                                        # Diference utilized just for simplicity
    bulk_perm = perm_init                                          # Initial permitivity = epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.01                                                      # Diferentiation step
                                                                   # Diferentiation until p = 1
    while x<1:                                                    
        dy = ((bulk_perm*diff)/(1-diff+x*diff)) * ((wp-bulk_perm)/(Lw*wp+(1-Lw)*bulk_perm))
        x=x+dx
        bulk_perm = bulk_perm+dy*dx
        
    return bulk_perm


def LR_MV(water, por, ap, sp, wp, CEC): 
    """
    Calculate the soil bulk real relative dielectric permittivity using the Lichtenecker and Rother model.

    This function computes the bulk real relative dielectric permittivity of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1], and Mendoza Veirana [2] alpha correction model (LR_MV). 
    The model incorporates the water content, bulk and particle densities, and 
    permittivities of air, solid, and water, as well as the soil cation exchange capacity. Reported R2 = 0.94

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    por: array_like
        Soil porosity [m**3/m**3].
    ap : array_like
        Soil air real relative dielectric permittivity phase [-].
    sp : array_like
        Soil solid real relative dielectric permittivity phase [-].
    wp : array_like
        Soil water phase real dielectric permittivity [-].
    CEC : array_like
        Soil cation exchange capacity [meq/100g].

    Returns
    -------
    array_like
        The estimated soil bulk real relative dielectric permittivity [-]

    References
    ----------
    .. [1] Lichtenecker, K., and K. Rother. 1931. 
    Die Herleitung des logarithmischen Mischungsgesetzes aus allgemeinen Prinzipien der staionären Strömung. 
    Phys. Z. 32:255-260.
    .. [2] Gaston Mendoza Veirana, Jeroen Verhegge, Wim Cornelis, Philippe De Smedt, 
    Soil dielectric permittivity modelling for 50-MHz instrumentation,
    Geoderma, Volume 438, 2023, 116624, ISSN 0016-7061,
    https://doi.org/10.1016/j.geoderma.2023.116624.
    
    Example
    -------
    >>> LR_MV(0.3, 1.3, 2.65, 1, 4, 80, 20)
    28.577  

    """
    alpha = 0.248*np.log(CEC) + 0.366
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def LR(water, por, ap, sp, wp, alpha): 
    """
    Calculate the soil bulk real relative dielectric permittivity using the Lichtenecker and Rother model.

    This function computes the bulk real relative dielectric permittivity of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1] (LR). 
    The model incorporates the water content, bulk and particle densities, and 
    permittivities of air, solid, and water.

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    por: array_like
        Soil porosity [m**3/m**3].
    ap : array_like
        Soil air real relative dielectric permittivity phase [-].
    sp : array_like
        Soil solid real relative dielectric permittivity phase [-].
    wp : array_like
        Soil water phase real dielectric permittivity [-].
    alpha : float
        Soil alpha exponent as defined in volumetric mixing theory [-].

    Returns
    -------
    array_like
        The estimated soil bulk real relative dielectric permittivity [-]

    References
    ----------
    .. [1] Lichtenecker, K., and K. Rother. 1931. 
    Die Herleitung des logarithmischen Mischungsgesetzes aus allgemeinen Prinzipien der staionären Strömung. 
    Phys. Z. 32:255-260.
    
    Example
    -------
    >>> LR(0.3, 1.3, 2.65, 1, 4, 80, 0.5)
    15.006

    """
    if not isinstance(alpha, np.floating):
        alpha = alpha[0]
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def LR_W(water, por, ap, sp, wp, clay): 
    """
    Calculate the soil bulk real relative dielectric permittivity using the Lichtenecker and Rother model.

    This function computes the bulk real relative dielectric permittivity of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1], and Wunderlich [2] alpha correction model (LR_W). 
    The model incorporates the water content, bulk and particle densities, and 
    permittivities of air, solid, and water, as well as the soil clay content.

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    por: array_like
        Soil porosity [m**3/m**3].
    ap : array_like
        Soil air real relative dielectric permittivity phase [-].
    sp : array_like
        Soil solid real relative dielectric permittivity phase [-].
    wp : array_like
        Soil water phase real dielectric permittivity [-].
    clay : array_like
        Soil clay content [g/g]*100

    Returns
    -------
    array_like
        The estimated soil bulk real relative dielectric permittivity [-]

    References
    ----------
    .. [1] Lichtenecker, K., and K. Rother. 1931. 
    Die Herleitung des logarithmischen Mischungsgesetzes aus allgemeinen Prinzipien der staionären Strömung. 
    Phys. Z. 32:255-260.
    .. [2] Tina Wunderlich, Hauke Petersen, Said Attia al Hagrey, Wolfgang Rabbel; 
    Pedophysical Models for Resistivity and Permittivity of Partially Water-Saturated Soils. 
    Vadose Zone Journal 2013;; 12 (4): vzj2013.01.0023. doi: https://doi.org/10.2136/vzj2013.01.0023
    
    Example
    -------
    >>> LR_W(0.3, 1.3, 2.65, 1, 4, 80, 20)
    17.505 

    """
    alpha = -0.46*(clay/100)+0.71
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def LongmireSmithP(bulk_ec_dc, bulk_perm_inf, frequency_perm):
    """
    Calculate the soil bulk real relative dielectric permittivity using the Longmire-Smith model.

    This is a semiempirical model that calculates the soil bulk real relative dielectric permittivity at different
    electromagnetic frequencies [1].

    Parameters
    ----------
    bulk_ec_dc : array_like
        Soil bulk real direct current electrical conductivity [-].
    bulk_perm_inf : array_like
        Soil bulk real relative permittivity at infinite frequency [-].
    frequency_perm : array_like
        Frequency of dielectric permittivity measurement [Hz].

    Returns
    -------
    array_like
        Soil bulk real relative dielectric permittivity [-].

    Notes
    -----
    The Longmire-Smith equation uses a set of coefficients to account for the 
    frequency-dependent dielectric dispersion. If all values in the `bulk_ec_dc` 
    array are zero, the function returns 'bulk_perm_inf'.

    Global Variables Used
    ---------------------
    epsilon_0 : float
        The vacuum permittivity constant.

    References
    ----------
    .. [1] K. S. Smith and C. L. Longmire, “A universal impedance for soils,” 
    Defense Nuclear Agency, Alexandria, VA, USA, Topical 
    Report for Period Jul. 1 1975-Sep. 30 1975, 1975.

    Example
    -------
    >>> LongmireSmithP(0.1, 5, 50e6)
    23.328

    """
    
    if (bulk_ec_dc == 0).all():
        return bulk_perm_inf
    
    a = [3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1, 4.8, 2.17, 9.8e-1, 3.92e-1, 1.73e-1]
    f = (125*bulk_ec_dc)**0.8312
    bulk_permi_ = []

    for i in range(len(a)):
        F_ = f*(10**i)
        bulk_permi = a[i]/(1+(frequency_perm/F_)**2)
        bulk_permi_.append(bulk_permi)
    bulk_perm = bulk_perm_inf + sum(bulk_permi_)

    return bulk_perm


def Hilhorst(bulk_ec, water_ec, water_perm, offset_perm):
    """
    Calculate the soil bulk real relative dielectric permittivity using the Hilhorst model and return

    This function calculates the soil bulk real relative dielectric permittivity of a 
    soil-water mixture based on Hilhorst's model. The relation 
    connects the bulk electrical conductivity of the mixture with the permittivity 
    of the water phase and an offset for the permittivity.

    Parameters
    ----------
    bulk_ec : array_like
        Soil bulk real relative dielectric permittivity [-].
    water_ec : array_like
        Soil water real electrical conductivity [S/m].
    water_perm : array_like
        Soil water phase real dielectric permittivity [-]. 
    offset_perm : array_like
        Soil bulk real relative dielectric permittivity when soil bulk real electrical conductivity is zero [-].

    Returns
    -------
    array_like
        Soil bulk real relative dielectric permittivity [-].

    References
    ----------
    .. [1] Hilhorst, M.A. (2000), A Pore Water Conductivity Sensor. 
    Soil Sci. Soc. Am. J., 64: 1922-1925. https://doi.org/10.2136/sssaj2000.6461922x   

    Example
    -------
    >>> Hilhorst(0.05, 0.5, 80, 4)
    12.0

    """
    bulk_perm = bulk_ec*water_perm/water_ec + offset_perm

    return bulk_perm