import numpy as np
from scipy.constants import pi, epsilon_0

def WunderlichEC(water, ec_init, wat_init, wc, Lw):  
    """
    Calculate the soil bulk real electrical conductivity using the Wunderlich model.

    This is a effective medium model that uses a differential 
    approach to compute the electrical conductivity based on the initial 
    conditions and the water content [1]. Reported RMSE = 0.0055

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    ec_init : float
        Initial soil bulk real electrical conductivity [S/m].
    wat_init : float
        Initial soil volumetric water content [m**3/m**3].
    wc : array_like
        Soil water real electrical conductivity [S/m].
    Lw : float
        Soil scalar depolarization factor of water aggregates (effective medium theory) [-]

    Returns
    -------
    array_like
        The estimated soil bulk real electrical conductivity [S/m].

    Notes
    -----
    The Wunderlich model is a differential model and the approach used here 
    employs a simple while loop with a step of 0.01 until the differential 
    fraction reaches 1. The method takes into account the initial water content, 
    initial electrical conductivity, and weighting factors to determine the 
    bulk electrical conductivity.

    References
    ----------
    .. [1] Tina Wunderlich, Hauke Petersen, Said Attia al Hagrey, Wolfgang Rabbel; 
    Pedophysical Models for Resistivity and Permittivity of Partially Water-Saturated Soils. 
    Vadose Zone Journal 2013;; 12 (4): vzj2013.01.0023. doi: https://doi.org/10.2136/vzj2013.01.0023

    Example
    -------
    >>> WunderlichEC(0.3, 0.08, 0.05, 0.50, 0.01)
    0.182634

    """
    diff = water - wat_init                                               # Diference utilized just for simplicity
    bulk_ec = ec_init                                                     # Initial permitivity = Epsilon sub 1  
    x = 0                                                                 # Diferentiation from p = 0  
    dx = 0.01                                                             # Diferentiation step

    while x<1:                                                            # Diferentiation until p = 1
        dy = dx*((bulk_ec*diff)/(1-diff+x*diff)) * ((wc-(bulk_ec))/(Lw*wc + (1-Lw)*bulk_ec))
        x=x+dx
        bulk_ec=bulk_ec+dy

    return bulk_ec    


def Fu(water, clay, bd, pd, wc, solid_ec, dry_ec, sat_ec, s=1, w=2):
    """
    Calculate the soil bulk real electrical conductivity using the Fu model.

    This is a volumetric mixing model that takes into account various soil properties 
    such as clay content, bulk density, particle density, and water content. 
    It was exhaustively validated using several soil samples [1]. Reported R2 = 0.98

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    clay : array_like
        Soil clay content [g/g]*100.
    bd : array_like 
        Soil bulk density (g/cm^3).
    pd : array_like
        Soil particle density (g/cm^3).
    wc : array_like
        Soil water real electrical conductivity [S/m].
    solid_ec : array_like
        Soil solid real electrical conductivity [S/m].
    dry_ec : array_like
        Soil bulk real electrical conductivity at zero water content [S/m].
    sat_ec : array_like
        Soil bulk real electrical conductivity at saturation water content [S/m].
    s : float, optional
        Phase exponent of the solid, default is 1.
    w : float, optional
        Phase exponent of the water, default is 2.

    Returns
    -------
    array_like
        The estimated bulk electrical conductivity [S/m].

    Notes
    -----
    The method uses default values for s and w, which are 1 and 2 respectively, 
    but can be modified if necessary. Three different forms of the model are used 
    depending on the soil data availability. The soil electrical conductivity of solid surfaces 
    is calculated as in [1] using the formula of Doussan and Ruy (2009) [2]

    References
    ----------
    .. [1] Yongwei Fu, Robert Horton, Tusheng Ren, J.L. Heitman,
    A general form of Archie's model for estimating bulk soil electrical conductivity,
    Journal of Hydrology, Volume 597, 2021, 126160, ISSN 0022-1694, https://doi.org/10.1016/j.jhydrol.2021.126160.
    .. [2] Doussan, C., and Ruy, S. (2009), 
    Prediction of unsaturated soil hydraulic conductivity with electrical conductivity, 
    Water Resour. Res., 45, W10408, doi:10.1029/2008WR007309.

    Example
    -------
    >>> Fu(0.3, 30, 1.3, 2.65, 0.3, 0, np.nan, np.nan)
    0.072626

    """
    d = 0.6539
    e = 0.0183
    por = 1 - bd/pd
    surf_ec = (d*clay/(100-clay))+e # Soil electrical conductivity of solid surfaces

    if np.isnan(dry_ec) & np.isnan(sat_ec):
        bulk_ec = solid_ec*(1-por)**s + (water**(w-1))*(por*surf_ec) + wc*water**w

    elif ~(np.isnan(dry_ec)) & ~(np.isnan(sat_ec)):
        bulk_ec = dry_ec + ((dry_ec-sat_ec)/(por**w) - surf_ec)*water**w + (water**(w-1))*(por*surf_ec)

    elif ~(np.isnan(dry_ec)) & np.isnan(sat_ec):
        sat_ec = dry_ec + (wc+surf_ec)*por**w
        bulk_ec = dry_ec + ((dry_ec-sat_ec)/(por**w) - surf_ec)*water**w + (water**(w-1))*(por*surf_ec)

    return bulk_ec


def LongmireSmithEC(bulk_ec_dc, frequency_ec):
    """
    Calculate the soil bulk real electrical conductivity using the Longmire-Smith model.

    This is a semiempirical model that calculates the soil bulk real electrical conductivity at different
    electromagnetic frequencies [1].

    Parameters
    ----------
    bulk_ec_dc : array_like
        Soil bulk real direct current electrical conductivity [S/m].
    frequency_ec : array_like
        Frequency of electric conductivity measurement [Hz].

    Returns
    -------
    array_like
        Soil bulk real electrical conductivity [S/m].

    Notes
    -----
    The Longmire-Smith equation uses a set of coefficients to account for the 
    frequency-dependent dielectric dispersion. If all values in the `bulk_ec_dc` 
    array are zero, the function returns 0.

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
    >>> LongmireSmithEC(np.array([0.05, 0.10]), 130)
    array([0.05153802, 0.10245936])

    """
    if (bulk_ec_dc == 0).all():
        return 0
    
    else: 
        a = [3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1, 4.8, 2.17, 9.8e-1, 3.92e-1, 1.73e-1]
        f = (125*bulk_ec_dc)**0.8312
        bulk_eci_ = []
        
        for i in range(len(a)):
            F_ = f*(10**i)
            bulk_eci = 2*pi*epsilon_0*(a[i]*F_*(frequency_ec/F_)**2/(1+(frequency_ec/F_)**2))                     
            bulk_eci_.append(bulk_eci)

        bulk_ec = bulk_ec_dc + sum(bulk_eci_)
        return bulk_ec
    
    
def Rhoades(water, wc, s_ec, E, F):
    """
    Calculate the bulk electrical conductivity using the Rhoades model.

    This function estimates the bulk electrical conductivity of a soil-water mixture 
    using the Rhoades equation [1]. The model combines the contributions of the water's 
    electrical conductivity and the soil's electrical conductivity, adjusted by 
    weighting factors based on the water content.

    Parameters
    ----------
    water : array_like
        Soil volumetric water content [m**3/m**3].
    wc : array_like
        Soil water real electrical conductivity [S/m].
    s_ec : array_like
        Soil bulk real surface electrical conductivity [S/m].
    E : float
        Weighting factor for the quadratic term of water content in the Rhoades equation.
    F : float
        Weighting factor for the linear term of water content in the Rhoades equation.

    Returns
    -------
    array_like
        Soil bulk real electrical conductivity [S/m].

    References
    ----------
    .. [1] Rhoades, J.D., Raats, P.A.C. and Prather, R.J. (1976), 
    Effects of Liquid-phase Electrical Conductivity, Water Content, 
    and Surface Conductivity on Bulk Soil Electrical Conductivity. 
    Soil Science Society of America Journal, 40: 651-655. 
    https://doi.org/10.2136/sssaj1976.03615995004000050017x

    Example
    -------
    >>> Rhoades(0.3, 0.5, 0.001, 1, 0.5)
    0.121

    """
    bulk_ec = wc*(E*(water**2)+F*water) + s_ec
    
    return bulk_ec