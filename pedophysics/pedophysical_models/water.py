import numpy as np

def LR_W(bp, bd, pdn, ap, sp, wp, clay): 
    """
        Lichtenecker and Rother, 1931 and Wunderlich et al., 2013
        
        Parameters
        ----------
        bulk_perm: float
            Real relative dielectric permittivity [-]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]
            
        ap: float
            air permittivity phase [-]
            
        sp: float
            solid permittivity phase [-]
            
        wp: float
            water permittivity phase [-]
            
        clay: float
            Soil volumetric clay content [m**3/m**3]
            
        Returns
        -------
        Volumetric water content: float
    """
    por = 1 - bd/pdn    
    alpha = -0.46*(clay/100)+0.71
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water


def LR(bp, bd, pdn, ap, sp, wp, alpha): 
    """
        Lichtenecker and Rother, 1931
        
        Parameters
        ----------
        bulk_perm: float
            Real relative dielectric permittivity [-]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]
            
        ap: float
            air permittivity phase [-]
            
        sp: float
            solid permittivity phase [-]
            
        wp: float
            water permittivity phase [-]

        alpha: float
            Soil alpha exponent as defined in volumetric mixing theory [-]
            
        Returns
        -------
        Volumetric water content: float
    """
    alpha = alpha [0]

    if np.isnan(alpha):
        alpha = 0.5
    por = 1 - bd/pdn    
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water


def LR_MV(bp, bd, pdn, ap, sp, wp, CEC): 
    """
        Lichtenecker and Rother, 1931 and Mendoza Veirana 
        
        Parameters
        ----------
        bulk_perm: float
            Real relative dielectric permittivity [-]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]
            
        ap: float
            air permittivity phase [-]
            
        sp: float
            solid permittivity phase [-]
            
        wp: float
            water permittivity phase [-]
            
        CEC: float
            Cation exchange capacity [meq/100g]

        Returns
        -------
        Volumetric water content: float
    """
    por = 1 - bd/pdn  
    alpha = 0.248*np.log(CEC) + 0.366  
    water = (bp**alpha - (1-por)*sp**alpha - por*ap**alpha) / (wp**alpha - ap**alpha)

    return water


#def fu(bulk_cond, clay, bd, pd, wec, sc, dry_c, sat_c, s=1):
    """
        Fu et al., 2021 
        
        Parameters
        ----------
        bulk_cond: float
            Soil bulk real electrical conductivity [S/m]

        clay: float
            Soil volumetric clay content [m**3/m**3]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]

        wec: float
            Soil water real electrical conductivity [S/m]
 
        sc: float
            Soil solid real electrical conductivity [S/m]

        w: float
            phase exponent of the water [-]

        s: float
            phase exponent of the solid [-]

        Returns
        -------
        water: float
            Volumetric water content [-]
    """      
    d = 0.6539
    e = 0.0183
    por = 1 - bd/pd
    surf_cond = d*clay/(100-clay)+e

    if np.isnan(dry_c) & np.isnan(sat_c):
        p = [wec, por*surf_cond, (sc*(1-por)**s)-bulk_cond]
        roots = np.roots(p)
        roots = roots[roots.imag == 0 ]
        roots = roots[roots > 0]
        water = roots.real
        print(water)

    elif ~(np.isnan(dry_c)) & ~(np.isnan(sat_c)):
        p = [(dry_c-sat_c)/(por**2) - surf_cond, por*surf_cond, dry_c-bulk_cond]
        roots = np.roots(p)
        roots = roots[roots.imag == 0 ]
        water = roots[roots > 0]
        water = roots.real

    elif ~(np.isnan(dry_c)) & np.isnan(sat_c):
        sat_c = dry_c + (wec+surf_cond)*por**2
        p = [(dry_c-sat_c)/(por**2) - surf_cond, por*surf_cond, dry_c-bulk_cond]
        roots = np.roots(p)
        roots = roots[roots.imag == 0 ]
        water = roots[roots > 0]
        water = roots.real

    return water