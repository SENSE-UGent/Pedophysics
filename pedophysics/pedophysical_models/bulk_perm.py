import numpy as np

def LongmireSmithP(bulk_cond, bulk_perm_inf, frequency_perm):
    """
        Longmire and Smith 1975 
        
        Parameters
        ----------
        bulk_cond: float
            Soil bulk real direct current electrical conductivity [S/m]
        
        bulk_perm_inf: float
            Bulk permittivity at infinite frequency [-]
            
        frequency_perm: float
            frequency of permittivity readings [-]
   
        Returns
        -------
        bulk_perm: float
            Soil bulk real relative dielectric permittivity       
    """ 
    a = [3.4e6, 2.74e5, 2.58e4, 3.38e3, 5.26e2, 1.33e2, 2.72e1, 1.25e1, 4.8, 2.17, 9.8e-1, 3.92e-1, 1.73e-1]
    f = (125*bulk_cond)**0.8312
    bulk_permi_ = []

    for i in range(len(a)):
        F_ = f*(10**i)
        bulk_permi = a[i]/(1+(frequency_perm/F_)**2)
        bulk_permi_.append(bulk_permi)

    bulk_perm = bulk_perm_inf + sum(bulk_permi_)
    return bulk_perm


def WunderlichP(water, perm_init, wat_init, wp, Lw): 
    """
        Wunderlich et.al 2013 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        perm_init: float
            minimum real permittivity [-]
            
        wat_init: float
            minimum volumetric water content [-]
            
        wp: float
            water permittivity phase [-]
            
        Lw: float
            Soil scalar depolarization factor of water agregates [-]
   
        Returns
        -------
        bulk_perm: float
            Soil bulk real relative dielectric permittivity   
    """     
    diff = water - wat_init                                          # Diference utilized just for simplicity
    bulk_perm = perm_init                                            # Initial permitivity = epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.01                                                      # Diferentiation step
                                                                   # Diferentiation until p = 1
    while x<1:                                                    
        dy = ((bulk_perm*diff)/(1-diff+x*diff)) * ((wp-bulk_perm)/(Lw*wp+(1-Lw)*bulk_perm))
        x=x+dx
        bulk_perm = bulk_perm+dy*dx
        
    return bulk_perm


def RothMV(water, bd, pdn, ap, sp, wp, CEC): 
    """
        Roth et al., 1990 and Mendoza Veirana et al., 2022
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
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
        bulk_perm: float
            Soil bulk real relative dielectric permittivity   
    """
    por = 1 - bd/pdn    

    alpha = 0.248*np.log(CEC) + 0.366
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def RothCRIM(water, bd, pdn, ap, sp, wp, alpha): 
    """
        Roth et al., 1990
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
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
            Alpha exponent as in Roth's model
            
        Returns
        -------
        bulk_perm: float
            Soil bulk real relative dielectric permittivity   
    """
    por = 1 - bd/pdn    
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def RothW(water, bd, pdn, ap, sp, wp, clay): 
    """
        Roth et al., 1990 and Wunderlich et al., 2013
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
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
        bulk_perm: float
            Soil bulk real relative dielectric permittivity   
    """

    por = 1 - bd/pdn    
    alpha = -0.46*(clay/100)+0.71
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm
