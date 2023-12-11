import numpy as np

def LR_W(bp, bd, pdn, ap, sp, wp, clay): 
    """
    Calculate the soil volumetric water content using the Lichtenecker and Rother model.

    This function computes the volumetric water content of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1], and Wunderlich [2] alpha correction model (LR_W). 
    The model incorporates the bulk dielectric permittivity, bulk and particle densities, and 
    permittivities of air, solid, and water, as well as the soil clay content.

    Parameters
    ----------
    bp : array_like
        Soil bulk real relative dielectric permittivity [-]
    bd : array_like 
        Soil bulk density (g/cm^3).
    pd : array_like
        Soil particle density (g/cm^3).
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
        Soil volumetric water content [m**3/m**3].

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
    >>> LR_W(10, 1.3, 2.65, 1, 4, 80, 20)
    0.177
    """

    por = 1 - bd/pdn    
    alpha = -0.46*(clay/100)+0.71
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water


def LR(bp, bd, pdn, ap, sp, wp, alpha): 
    """
    Calculate the soil volumetric water content using the Lichtenecker and Rother model.

    This function computes the volumetric water content of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1] (LR). 
    The model incorporates the bulk dielectric permittivity, bulk and particle densities, and 
    permittivities of air, solid, and water.

    Parameters
    ----------
    bp : array_like
        Soil bulk real relative dielectric permittivity [-]
    bd : array_like 
        Soil bulk density (g/cm^3).
    pdn : array_like
        Soil particle density (g/cm^3).
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
        Soil volumetric water content [m**3/m**3].

    References
    ----------
    .. [1] Lichtenecker, K., and K. Rother. 1931. 
    Die Herleitung des logarithmischen Mischungsgesetzes aus allgemeinen Prinzipien der staionären Strömung. 
    Phys. Z. 32:255-260.
    
    Example
    -------
    >>> LR(10, 1.3, 2.65, 1, 4, 80, 0.5)
    0.210

    """
    if not isinstance(alpha, np.floating):
        alpha = alpha[0]

    if np.isnan(alpha):
        alpha = 0.5
    por = 1 - bd/pdn    
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water


def LR_MV(bp, bd, pdn, ap, sp, wp, CEC): 
    """
    Calculate the soil volumetric water content using the Lichtenecker and Rother model.

    This function computes the volumetric water content of a soil 
    mixture using the volumetric mixing model of Lichtenecker and Rother [1], and Mendoza Veirana [2] alpha correction model (LR_MV). 
    The model incorporates the bulk dielectric permittivity, bulk and particle densities, and 
    permittivities of air, solid, and water, as well as the soil cation exchange capacity. Reported R2 = 0.94

    Parameters
    ----------
    bp : array_like
        Soil bulk real relative dielectric permittivity [-]
    bd : array_like 
        Soil bulk density (g/cm^3).
    pdn : array_like
        Soil particle density (g/cm^3).
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
        Soil volumetric water content [m**3/m**3].

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
    >>> LR_MV(10, 1.3, 2.65, 1, 4, 80, 20)
    0.078

    """
    por = 1 - bd/pdn  
    alpha = 0.248*np.log(CEC) + 0.366  
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water