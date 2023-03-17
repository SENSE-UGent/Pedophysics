import numpy as np
from scipy.constants import pi, epsilon_0

def WunderlichC(water, cond_init, wat_init, wc, Lw):  
    """
        Wunderlich et.al 2013 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        cond_init: float
            minimum electrical conductivity [S/m]
            
        wat_init: float
            minimum volumetric water content [-]
            
        wc: float
            Soil water real electrical conductivity [-]
            
        Lw: float
            Soil scalar depolarization factor of water agregates [-]
   
        Returns
        -------
        bulk_cond: float
            Soil bulk real electrical conductivity [S/m]
    """     
    diff = water - wat_init                                               # Diference utilized just for simplicity
    bulk_cond = cond_init                                                 # Initial permitivity = Epsilon sub 1  
    x = 0                                                                 # Diferentiation from p = 0  
    dx = 0.01                                                            # Diferentiation step

    while x<1:                                                            # Diferentiation until p = 1
        dy = dx*((bulk_cond*diff)/(1-diff+x*diff)) * ((wc-(bulk_cond))/(Lw*wc + (1-Lw)*bulk_cond))
        x=x+dx
        bulk_cond=bulk_cond+dy

    return bulk_cond    


def Fu(water, clay, bd, pd, wc, sc, dry_c, sat_c, s=1, w=2):
    """
        Fu et al., 2021 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]

        clay: float
            Soil volumetric clay content [m**3/m**3]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]

        wc: float
            Soil water real electrical conductivity [S/m]
 
        sc: float
            Soil solid real electrical conductivity [S/m]

        w: float
            phase exponent of the water [-]

        s: float
            phase exponent of the solid [-]

        Returns
        -------
        bulk_cond: float
            Soil bulk real electrical conductivity [S/m]
    """      
    d = 0.6539
    e = 0.0183
    por = 1 - bd/pd
    surf_cond = d*clay/(100-clay)+e

    if np.isnan(dry_c) & np.isnan(sat_c):
        bulk_cond = sc*(1-por)**s + wc*water**w + (water**(w-1))*(por*surf_cond)

    elif ~(np.isnan(dry_c)) & ~(np.isnan(sat_c)):
        bulk_cond = dry_c + ((dry_c-sat_c)/(por**w) - surf_cond)*water**w + (water**(w-1))*(por*surf_cond)

    elif ~(np.isnan(dry_c)) & np.isnan(sat_c):
        sat_c = dry_c + (wc+surf_cond)*por**w
        bulk_cond = dry_c + ((dry_c-sat_c)/(por**w) - surf_cond)*water**w + (water**(w-1))*(por*surf_cond)

    return bulk_cond


def LongmireSmithC(bulk_cond, frequency_cond):
    """
        Longmire and Smith 1975 
        
        Parameters
        ----------
        bulk_cond: float
            Soil bulk real direct current electrical conductivity [-]
            
        frequency_cond: float
            frequency of electrical conductivity readings [-]
   
        Returns
        -------
        bulk_perm: float
            Soil bulk real relative dielectric permittivity       
    """               
    if bulk_cond == 0:
        return 0
    
    a = [3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1, 4.8, 2.17, 9.8e-1, 3.92e-1, 1.73e-1]
    f = (125*bulk_cond)**0.8312
    bulk_condi_ = []
    
    for i in range(len(a)):
        F_ = f*(10**i)
        bulk_condi = 2*pi*epsilon_0*(a[i]*F_*(frequency_cond/F_)**2/(1+(frequency_cond/F_)**2))                     
        bulk_condi_.append(bulk_condi)

    bulk_cond_f = bulk_cond + sum(bulk_condi_)
    return bulk_cond_f