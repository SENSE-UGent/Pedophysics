import numpy as np
from scipy.constants import pi, epsilon_0

def WunderlichEC(water, ec_init, wat_init, wc, Lw):  
    """
        Wunderlich et.al 2013 
        
        Parameters
        ----------
        water: float
            volumetric water content [-]
        
        ec_init: float
            minimum electrical conductivity [S/m]
            
        wat_init: float
            minimum volumetric water content [-]
            
        wc: float
            Soil water real electrical conductivity [-]
            
        Lw: float
            Soil scalar depolarization factor of water agregates [-]

        Returns
        -------
        bulk_ec: float
            Soil bulk real electrical conductivity [S/m]
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
        Fu et al., 2021 
        
        Parameters
        ----------
        water: float
            volumetric water content [-]

        clay: float
            Soil volumetric clay content [m**3/m**3]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]

        wc: float
            Soil water real electrical conductivity [S/m]
 
        solid_ec: float
            Soil solid real electrical conductivity [S/m]

        dry_ec: float
            Soil bulk real electrical conductivity with zero water content [S/m]

        sat_ec: float
            Soil bulk real electrical conductivity at saturation capacity [S/m]

        w: float
            phase exponent of the water [-]

        s: float
            phase exponent of the solid [-]

        Returns
        -------
        bulk_ec: float
            Soil bulk real electrical conductivity [S/m]
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


def LongmireSmithEC(bulk_ec, frequency_ec):
    """
        Longmire and Smith, 1975 
        
        Parameters
        ----------
        bulk_ec: float
            Soil bulk real direct current electrical conductivity [-]
            
        frequency_ec: float
            frequency of electrical conductivity readings [-]
   
        Returns
        -------
        bulk_perm: float
            Soil bulk real relative dielectric permittivity       
    """               
    if (bulk_ec == 0).all():
        return 0
    
    else: 
        a = [3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1, 4.8, 2.17, 9.8e-1, 3.92e-1, 1.73e-1]
        f = (125*bulk_ec)**0.8312
        bulk_eci_ = []
        
        for i in range(len(a)):
            F_ = f*(10**i)
            bulk_eci = 2*pi*epsilon_0*(a[i]*F_*(frequency_ec/F_)**2/(1+(frequency_ec/F_)**2))                     
            bulk_eci_.append(bulk_eci)

        bulk_ec_f = bulk_ec + sum(bulk_eci_)
        return bulk_ec_f
    
    
def Rhoades(water, water_ec, s_ec, E, F):
    """
        Rhoades et al., 1976
        
        Parameters
        ----------
        water: float
            volumetric water content [-]

        water_ec: float
            Soil water real electrical conductivity [S/m]
 
        s_ec: float
            Soil bulk real surface electrical conductivity [S/m]

        E: float
            Empirical constant [-]

        f: float
            Empirical constant [-]

        Returns
        -------
        bulk_ec: float
            Soil bulk real electrical conductivity [S/m]
    """  
    bulk_ec = water_ec*(E*(water**2)+F*water) + s_ec
    
    return bulk_ec