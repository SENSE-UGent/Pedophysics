"""
    Pedophysics
    ...

    :AUTHOR: Gaston Mendoza Veirana
    :CONTACT: gaston.mendozaveirana@ugent.be

    :REQUIRES: numpy, scipy
"""

# Import
import numpy as np
from scipy.constants import epsilon_0, pi

def longmire_smith(bulk_cond, bulk_perm_inf, frequency_perm):
    """
        Longmire and Smith 1975 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        bulk_perm_inf: float
            Bulk permittivity at infinite frequency [-]
            
        frequency_perm: float
            frequency of permittivity readings [-]
   
        Returns
        -------
        Real bulk permittivity: float   
    """               
    a1=3.4e6
    a2=2.74e5
    a3=2.58e4
    a4=3.38e3
    a5=5.26e2
    a6=1.33e2
    a7=2.72e1
    a8=1.25e1
    a9=4.8
    a10=2.17
    a11=9.8e-1
    a12=3.92e-1
    a13=1.73e-1
    f = (125*bulk_cond)**0.8312
    bulk_perm = bulk_perm_inf + a1/(1 + (frequency_perm/(f*(10**0)))) + a2/(1 + (frequency_perm/(f*(10**1)))) + a3/(1+ (frequency_perm/(f*(10**2)))) + a4/(1+ (frequency_perm/(f*(10**3)))) + a5/(1 + (frequency_perm/(f*(10**4)))) + a6/(1 + (frequency_perm/(f*(10**5)))) + a7/(1 + (frequency_perm/(f*(10**6)))) + a8/(1 + (frequency_perm/(f*(10**7)))) + a9/(1 + (frequency_perm/(f*(10**8)))) + a10/(1 + (frequency_perm/(f*(10**9)))) + a11/(1 + (frequency_perm/(f*(10**10)))) + a12/(1 + (frequency_perm/(f*(10**11)))) + a13/(1 + (frequency_perm/(f*(10**12))))
    
    return bulk_perm


def wunderlich_p(water, perm_init, wat_init, wp, L): 
    """
        Wunderlich et.al 2013 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        perm_init: float
            lowest real permittivity [-]
            
        wat_init: float
            lowest volumetric water content [-]
            
        wp: float
            water permittivity phase [-]
            
        L: float
            depolarization factor [-]
   
        Returns
        -------
        Real bulk permittivity: float   
    """     
    diff = water - wat_init                                          # Diference utilized just for simplicity
                                                                   # Initializing diferenciation parameters
    bulk_perm = perm_init                                                  # Initial permitivity = epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.01                                                      # Diferentiation step
                                                                   # Diferentiation until p = 1
    while x<1:                                                    
        dy = ((bulk_perm*diff)/(1-diff+x*diff)) * ((wp-bulk_perm)/(L*wp+(1-L)*bulk_perm))
        x=x+dx
        bulk_perm = bulk_perm+dy*dx
        
    return bulk_perm


def roth_mv(water, bd, pdn, ap, sp, wp, CEC): 
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
        Real bulk permittivity: float
    """
    por = 1 - bd/pdn    
    alpha = 0.248*np.log(CEC) + 0.366
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


def roth_wn(water, bd, pdn, ap, sp, wp, clay): 
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
        Real bulk permittivity: float
    """
    por = 1 - bd/pdn    
    alpha = -0.46*(clay/100)+0.71
    bulk_perm = ( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha)

    return bulk_perm


################################## water_perm ##############################################


def Malmberg_Maryott56(T):
    """
        Malmberg & Maryott56, 1956
        
        Parameters
        ----------
        temperature: float
            Soil temperature in Kelvin [K]
 
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]: float   
    """    
    Tc = T - 273.15 # Kelvin to Celsius
    water_perm = 87.740 - 0.40008*Tc + 9.398e-4*Tc**2 - 1.410e-6*Tc**3

    return water_perm


def olhoeft(T, c):
    """
        Olhoeft, 1986
        
        Parameters
        ----------
        temperature: float
            Soil temperature in Kelvin [K]
        
        salinity: float
            Soil (NaCl) salinity of the bulk pore fluid [mol/L]
 
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]: float   
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


def stogryn(T, c, f):
    """
        Stogryn, 1971
        
        Parameters
        ----------
        temperature: float
            Soil temperature in Kelvin [K]
        
        salinity: float
            Soil (NaCl) salinity of the bulk pore fluid [mol/L]
 
        frequency_perm: float
            frequency of permittivity readings [-]
            
        Returns
        -------
        water_perm: float
            Soil water phase real dielectric permittivity [-]: float      
    """
    molarity = 58.44 # Molar mass of NaCl
    water_perm_inf = 4.5

    c_ppt = c/molarity
    N = c_ppt*(1.707e-2 + 1.205e-5*c_ppt + 4.058e-9*c_ppt**2)
    a_N = 1 - 0.2551*N + 5.151e-2*N**2 - 6.889e-3*N**3
    b_N_T = 0.1463e-2*N*T + 1 - 0.04896*N - 0.02967*N**2 + 5.644e-3*N**3
    e_t_0 = 87.74 + 4.0008*T + 9.398e-4*T**2 + 1.41e-6*T**3
    two_pi_tau_T_0 = 1.1109e-10 - 3.824e-12*T + 6.938e-14*T**2 - 5.096e-16*T**3
    water_perm = water_perm_inf + (e_t_0*a_N - water_perm_inf)/(1 + (two_pi_tau_T_0*b_N_T*f)**2)

    return water_perm


################################ ELECTRICAL CONDUCTIVITY EQUATIONS #######################
#################################### UNCHECKED ###########################################
    

def parallel(water, bd, pdn, cond_matrix, pore_ec, cond_air = 0.01e-3):
    """
    Parallel model, Glover 2015 Table 2
    """
    por = 1 - bd/pdn
    
    
    cond = (1 - por)*cond_matrix + water*pore_ec+ (por - water)*cond_air
    return cond



def perpendicular(water, bd, pdn, cond_matrix, pore_ec, cond_air = 0.01e-3):
    """
    Perpendicular model, Glover 2015 Table 2
    """
    por = 1 - bd/pdn
    
    cond = (((1-por)/cond_matrix) + water/pore_ec+ (por - water)/cond_air)**-1
    return cond



def random(water, bd, pdn, cond_matrix, pore_ec, cond_air = 0.01e-3):
    """
    Random model, Glover 2015 Table 2
    """
    por = 1 - bd/pdn
    
    cond = (cond_matrix**(1-por)) * (pore_ec**water) * (cond_air**(por-water))
    return cond


def HSme(water, bd, pdn, cond_matrix, pore_ec):
    """Hashin and Shtrikman (1962) Commonly denoted as HS-. Derived from effective medium considerations. Glover 2015 Table2"""
    por = 1 - bd/pdn
    
    return cond_matrix*(1+(3*por*((pore_ec)-cond_matrix)/(3*cond_matrix+(1-por)*((pore_ec)-cond_matrix))))
       
                       
def HSma(water, bd, pdn, cond_matrix, pore_ec):
    """Hashin and Shtrikman (1962) Commonly denoted HS+. Derived from effective medium considerations. Glover 2015 Table2"""
    por = 1 - bd/pdn
    
    return (pore_ec)*(1-(3*(1-por)*((pore_ec)-cond_matrix)/((3/pore_ec)-(por)*((pore_ec)-cond_matrix))))


def LRmodel(water, bd, pdn, cond_matrix, pore_ec, m):    
    """Lichtenecker and Rother (1931), Korvin (1982) Derived from the theory of functional equations under appropriate
boundary conditions. Formally, the same as Archie’s law if s1=0. Glover 2015 Table 2"""
    por = 1 - bd/pdn
    
    return (((pore_ec)**(1/m))*water + (cond_matrix**(1/m))*(1-por) )**m


def archie(bd, pdn, pore_ec, m):
    """
    Achie original formulation for saturated sandstones
    m: cementation exponent
    """
    por = 1 - (bd/pdn)
    
    return ((por**-m)*pore_ec)**-1


def comb_archie(water, bd, pdn, pore_ec, m, n):
    """ Archie first and second law combination"""
    por = 1 - bd/pdn
    S = water/por                                        # Saturation 
   
    return (por**(m))*(S**(n))*(pore_ec)
                       

def glover10(water, bd, pdn, cond_matrix, pore_ec, m):
    """ Glover (2010b) Derived from the conventional Archie’s law by considering boundary conditions implied by geometric constraints. Tabla2, Glover 2015. Three phases"""
    por = 1 - bd/pdn
    p = np.log(1-por**m)/np.log(1-por)
    #m2 = (water**m)/water
    m2 = np.log(1-water**m)/np.log(1-water)   # This exponent needs also to consider the exponent of the air phase (neglected)
    
    return (pore_ec)*water**m + cond_matrix*(1-por)**m2


def shahsingh(water, clay_cont, pore_ec, c=1.45, m=1.25):
    """Shah and D. Singh, "Generalized Archie's Law for Estimation of Soil Electrical Conductivity"""
    if (clay_cont >= 5):
        c = .6 * clay_cont ** .55
        m = .92 * clay_cont ** .2
        
    return (c) * (pore_ec) * water ** m


def saey09H(clay_cont, organic_cont, LT):
    """ 
    Saey et al. 2009 (Ghent University) conductivity for horizontal coil horientation.
    S is set as unsaturated.
    Take a look to the original paper for details about how was measured clay_cont.
    
    LT = land use. This variable is dependent of organic_content, terefore must be set as = 1 if landuse = arable OR 
    = 0 if landuse = pasture. 
    """
    return (7.2 + 1.098 * clay_cont + 1.494 * LT * organic_cont + 1.339 * (1 - LT) * organic_cont)


def saey09V(clay_cont):
    """ 
    Saey et al. 2009 (Ghent University) conductivity for vertical coil horientation.
    S is set as unsaturated.
    Take a look to the original paper for details about how was measured clay_cont.
    """
    
    return (16.2 + 1.314 * clay_cont/100)


def peplinskicond(water, bd, pdn, clay_cont, sand_cont):
    """
    Empirical direct current electrical conductivity in imaginary component of complex permittivity of Peplinski 1995.  
    """
    alpha = 0.65
    
    beta = 1.33797 - 0.603*sand_cont/100 - 0.166*clay_cont/100
    condeff = 0.0467 + 0.2204*bd - 0.4111*sand_cont/100 + 0.6614*clay_cont/100
    
    cond = ((condeff)*(pdn - bd)*water**(beta/alpha))/(pdn*water)
    return cond


def rhoades90(water, bd, pdn, clay_cont, pore_ec):
    """
    Rhoades 1990 et al. 1b equation
    """
    thetas  =    bd/pdn
    thetaws =    0.639*water + 0.011 
    cond_solid = (0.023*clay_cont      - 0.021)
    
    return (((thetas + thetaws)**2)/thetas)*cond_solid + (water - thetaws)*pore_ec


def waxsmidt68(water, bd, pdn, CEC, pore_ec, m, B0):
    """Waxman & Smidt 1968"""
    e = 2.718280828
    B = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por = 1 - bd/pdn
    F = por**-m
    S = water/por   

    Qv = pdn*((1-por)/por)*CEC  
    surfc = B*Qv*S/F

    return (S**2)*pore_ec/F   +  surfc


def linde06(water, bd, pdn, CEC, pore_ec, m, n, B0):
    """
    Linde et al. 2006 conductivity equation. Surface conductivity proportion calculated as B*Qv    """
    e = 2.718280828
    B = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por = 1 - bd/pdn
    S = water/por                                        # Saturation 
    Qv = pdn*((1-por)/por)*CEC 
    cond_sup = B*Qv
    
    return (por**m) * ((S**n)* (pore_ec) + ((por**-m) - 1) * cond_sup)


def revil98(water, bd, pdn, CEC, pore_ec, m, n, B0, a=.4):
    """
    Revil et al 1998    """
    
    B = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por = 1 - bd/pdn
    S = water/por                                        # Saturation 
    Qv = pdn*((1-por)/por)*CEC                         # Excess of surface charge per unit volume
    
    return ((por**m)*(S**n)/a) * ((pore_ec)+(B*Qv/S))


def bardonpied(water, bd, pdn, clay_cont, cond_clay, pore_ec, m):
    """
    Bardon and Pied 1969 in Table 4 Glover 2015"""
    
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    F = por**-m
    
    return ((S**2)*pore_ec/F) + clay_cont*cond_clay*S


def schl(water, bd, pdn, clay_cont, cond_clay, pore_ec, m):
    """ 
    schlumberger eq in Table 4 Glover 2015 """
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    F = por**-m
    
    return (S**2)*pore_ec/(F*(1-clay_cont)) + clay_cont*cond_clay*S


def hossin(water, bd, pdn, clay_cont, cond_clay, pore_ec, m):
    """ Hossin 1960 in Table 4 Glover 2015"""
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    F = por**-m
    
    return ((S**2)*pore_ec/F) + (clay_cont**2)*cond_clay


def juhasz(water, bd, pdn, clay_cont, cond_clay, pore_ec, m, pdn_clay=2.85):
    """ 
    Juhasz 1981 in Table 4 Glover 2015"""
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    F = por**-m
    
    por_clay = 1 - (bd-0.2)/pdn_clay
    m_clay = np.log(cond_clay*pore_ec)/np.log(por_clay)
    F_clay = por_clay**-m_clay
        
    return ((S**2)*pore_ec/(F)) + ((cond_clay/F_clay) - pore_ec)*((clay_cont*por_clay*S)/(por)) 
    
    
def wund13c(water, L, pore_ec, cond_init, wat_init):  
    """
    #Wunderlich et.al 2013 Equation 22. This is trained fitting L, rest wat, wat_init and rest_init?    """
    dif = water - wat_init                                              # Diference utilized just for simplicity
    y = cond_init                                                     # Initial permitivity = Epsilon sub 1  
    x = 0                                                             # Diferentiation from p = 0  
    dx = 0.1                                                          # Diferentiation step
    
    while x<1:                                                        # Diferentiation until p = 1
        dy = dx*((((y)*(dif))/(1-dif+x*dif)) * (((pore_ec)-(y))/((L/pore_ec) + (1-L)*y)))
        x=x+dx
        y=y+dy
        
    return y                                                          # Return electrical conductivity of the soil


def logseq7(cond, freq, freqc, n):
    """Equation 7 in Logsdon 2010 et al. for real electrical conductivity vs frequency"""
    
    return cond*(1 + freq/freqc)**n


def logs10MHz(water, cond, freq):
    """Equation 7 in Logsdon 2010 et al. for real electrical conductivity vs frequency
       'n' exponent and central frequency is calculated using Table 2 Logsdon 2005"""
    freqc = 8.53  - 25.4 *water + 0.782*cond*1000
    n =     0.682 + 0.371*water
    
    return  cond*(1 + ((freq)/(freqc*1000000)))**n
    
    
def clavier77(water, bd_eff, pd_eff, clay_cont, clay_cond, clay_bd, pore_ec, m, n, clay_pd = 2.8, a=1):
    """Clavier et al. (1977) dual-water model as defined in Glover etal. 2015 pag 20, eq. 45"""
    por_tsh = 1 - (clay_bd/clay_pd)
    por_eff = (1 - (bd_eff/pd_eff))*(1-clay_cont/100) + por_tsh*clay_cont/100
    por_t =   por_eff + (clay_cont/100)*por_tsh
    S_bw =    (clay_cont/100)*por_tsh/por_t                    # Bound water saturation
    ro_bw =   (por_tsh**2)/clay_cond                           # Bound water conductivity
    S_wt =    S_bw + (water/por_eff)
                                                
    return    ((S_wt**n) - S_wt*(S_bw*(1-(a*ro_bw/pore_ec)))) / (a*pore_ec*por_t**m)
    
    
def mcbratney05(water, bd, pd, clay, cec):
    """ Semiempirical model of McBratney et al. 2005 adjusted to Edgeroi soil database (n> 1900) for APARENT electrical conductivity.
    cec is in mmol/g"""
    k = 631.8   # mS/m 
    por = 1 - bd/pd
    
    eca = k*(clay/100)*(water/por)*(cec*10)
    return eca/1000


def doussan_ruy(water, bd, pdn, sand, silt, clay, pore_ec, m):
    """Doussan & Ruy 2009. Surface conductivity based on WS model obtained for clay contents higher than 6 %. Based on Rhoades model otherwise. As defined in Doussan & Ruy 2009"""

    por = 1 - bd/pdn
    F = por**-m
    S = water/por   

    
    if clay > 6:
        scond = (0.654*clay/(sand+silt)) + 0.018
    else:
        scond = 0.023*clay - 0.0209

        
    return ((S**2)/F)*(pore_ec + scond/S)
   
    
def fu21(water, bd, pd, sand, silt, clay, cond_sat, cond_dry, w=2):
    """Fu et al. 2021"""
    
    por = 1 - bd/pd
    
    if clay > 6:
        scond = (0.654*clay/(sand+silt)) + 0.018
    else:
        scond = 0.023*clay - 0.0209
    
    return cond_dry + (((cond_sat-cond_dry)/por**w)-scond)*water**w + scond*por*water**(w-1)
    
    
    
################################### ELECTRICAL CONDUCTIVITY VS TEMPERATURE ###################################


def mcneill80(t, cond25, beta=0.02):
    """
    Mcneill 1981 equation
    """
    return (cond25*(1+beta*(t-25)))/1000


def arps53(t, cond25):
    """
    Arps 1953 Equation
    """
    return (cond25*((t+21.5)/46.5))/1000


def arps53_cond25(t, cond):
    """
    Arps 1953 Equation
    """
    return cond*(46.5/(t+21.5))


def luck09(t, cond25):
    """
    Luck 2009 Equation
    """
    return (cond25/(0.36+ np.exp((12.5-t)/28.5)))/1000


################################### ELECTRICAL CONDUCTIVITY OF WATER ###################################


def hilhorst_ecw(bulkcond, bulkperm, wp, offset):
    """ Hilhorst 2000 formula for water conductivity. It is given in the same utities as condbulk"""
    return (wp/bulkperm)*bulkcond + offset


def kelleners(water, bd, pdn, waterdens, condsat):
    """ Unknown deduction process. It is given in the same utities as condsat"""
    por = 1 - bd/pdn
    
    return (bd/waterdens)*(water/por)*(condsat/water)


def vogeler96(water, bulkcond):
    """ EC-TDR based measuremets, Water conductivity was determined by extraction process """
    return (bulkcond - (0.228*water - 0.042))/(0.804*water-0.217)
    
    
def weast_KCl(C):
    """"For a KC1 solution with concentrations ranging from 0.005 to 0.05 M, Weast (1965, p. D-81"""
    return (C+0.7e-4)/(7.6e-2)


def weast_KCl(C):
    """"For a CaC12 solution with concentrations ranging from 0.005 to 0.05 M, Weast (1965, p. D-81"""
    return (C+1.2e-3)/(9.4e-2)
       
    
def kargas17(bulkcond, ap_perm, wp, offset=6):
    """Kargas et al 2017 (Prediction of Soil Solution Electrical Conductivity by the Permittivity Corrected Linear Model Using a Dielectric Sensor) using Robinson et al 99 and Hilhorst formula """
    return (wp/(((ap_perm**0.5) - 0.628*bulkcond)**2) - offset)*bulkcond
    
    
def leao(ap_perm, ap_cond, sand):
    """"New semi-empirical formulae for predicting soil solution conductivityfrom dielectric properties at 50 MHzTairone P. Leao et al. 2010"""
    sand = sand * 100
    
    return (ap_cond-0.08)/((ap_perm-6.2)*(0.0057+0.000071*sand))


def bouksila(ap_perm, ap_cond, wp):
    """Soil water content and salinity determination using different dielectric methods in saline gypsiferous soil / Détermination de la teneur en eau et de la salinité de sols salins gypseux à l'aide de différentes méthodes diélectriques FETHI BOUKSILA , MAGNUS PERSSON , RONNY BERNDTSSON & AKISSA BAHRI"""
    
    offset = 0.4414*(ap_cond/10)**3 - 4.3435*(ap_cond/10)**2 + 13.733*ap_cond/10 + 3.9181
    
    return wp*ap_cond/(ap_perm - offset)
    
    
def fu21_cond_water(bd, pd, sand, silt, clay, cond_sat, cond_dry, w=2):
    """Fu et al. 2021"""
    
    por = 1 - bd/pd
    
    if clay > 6:
        scond = (0.654*clay/(sand+silt)) + 0.018
    else:
        scond = 0.023*clay - 0.0209
    
    return ((cond_sat - cond_dry)/por**w) - scond
    
    
############################   MODIFIED MODELS FOR ELECTRICAL CONDUCTIVITY    ######################################################################
    
    
def modcomb_archie(water, bd, pdn, clay, pore_ec, ps, n, offset):
    """ Archie first and second law combination. m calculed using prop 12"""
    
    alpha = (-0.46*clay/100)+0.71
    por = 1 - bd/pdn
    m = (np.log((((1 - por)*(ps/80)**alpha) + por)**(1/alpha))-offset/80) / np.log(por)
    S = water/por                                      # Saturation 
    
    return (por**(m))*(S**(n))*(pore_ec,)


def mod2comb_archie(water, bd, pdn, clay, pore_ec, ps, wp, offset):
    """ Archie first and second law combination. m calculed using prop 12"""
    
    alpha = (-0.46*clay/100)+0.71
    por = 1 - bd/pdn
    
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha))-offset/wp) / np.log(por)
    S = water/por                                      # Saturation 
    n = m
    
    return (por**(m))*(S**(n))*(pore_ec,)


def toppmod1(water, pore_ec,wp):
    "Modified Topp equation through Brovelli and Cassiani 2008 procedure. This did not work well"
    p = [4.3e-6*(wp**3)*(pore_ec**3), 
         -5.5e-4*(wp**2)*(pore_ec**2), 
         2.92e-2*(wp)*(pore_ec),
         -5.3e-2 - water]
    
    roots = np.roots(p) 
    roots = roots[roots.imag == 0 ]
    cond = roots[roots > 0]
    return cond[0].real
    
    
def toppmod2(water, pore_ec, wp, offset = 4):
    "Modified Topp equation through Hillhost 2000 formula. This did not work well"
    a = -5.3e-2
    b = 2.92e-2 
    c = -5.5e-4
    d = 4.3e-6
    
    p = [d*(wp**3)*(pore_ec**3), 
         2*d*(wp**2)*(pore_ec**2) + d*offset*(wp**2)*(pore_ec**2) + c*(wp**2)*(pore_ec**2), 
         2*d*(wp)*(pore_ec)*offset + d*(wp)*(pore_ec)*offset**2 + b*(wp)*(pore_ec) + 2*c*offset*(wp)*(pore_ec),
         d*offset**3 + offset**2 + b*offset + a - water]
    
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    cond = roots[roots > 0]
    return cond[0].real


def modlinde06(water, bd, pd, clay_cont, CEC, pore_ec,B0):
    """Linde et al. 2006 conductivity equation. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below), the effect seems negligible. 
    """
    if [clay_cont >= 5]:
        m = 0.92*clay_cont**0.2
        
    else:
        m = 1.25
        
    n = m    
    por = 1 - bd/pd
    F =    por**-m
    e    = 2.718280828
    
    B    = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    Qv   = pd*((1-por)/por)*CEC 
    cond_sup = B*Qv/F
    S = water/(por)
    
    return  (por**m) * (S**n)* (pore_ec,) + (1-(por**m))*B*Qv


def mod2linde06(water, bd, clay_cont, org, pore_ec,air_perm, ps, wp, n, cond_sup = .7, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """Linde et al. 2006 conductivity equation. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below), the effect seems negligible. 'm' is calculated using prop13 with a fixed offset==3. Surface conductivity is also fixed as 0.7 m S/m"""
    clay     = clay_cont/100
    org      = org/100
    somr     = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd       = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1

    por      = 1 - bd/pd
    alpha    = (-0.46 * clay) + 0.71
    S        = water/(por)
    offset   = 3                                   
    rate     = ps / wp
    
    #m        = np.log(((((1-por)*ratnp.expalpha) + water + ((por-water)*(1/wp)**alpha))**(1/alpha)- (offset/wp))*S**(-n))/ np.log(por)
    m= 1.50

    return     (por**m) * ((S**n)* (pore_ec,) + ((por**-m) - 1) * cond_sup/1000)


def mod3linde06(water, bd, clay_cont, org, pore_ec,air_perm, ps, wp, n, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """     Linde et al. 2006 conductivity equation. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below), the effect seems negligible. 'm' is calculated using prop13 with a fixed offset==4. Surface
     conductivity is also fixed as 1 m S/m 
     """
    clay     = clay_cont/100
    org      = org/100
    somr     = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd       = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    
    por      = 1 - bd/pd
    alpha    = (-0.46 * clay) + 0.71
    S        = water/(por)
    offset   = 4
    rate     = ps / wp
    
    #m        = np.log(((((1-por)*ratnp.expalpha) + water + ((por-water)*(1/wp)**alpha))**(1/alpha)- (offset/wp))*S**(-n))/ np.log(por)
    #print("Linde 3 m=", m)
    m = 1.5
    cond_sup = 1/1000
    #print("Mod 3 Linde cond sup=", cond_sup)
    return     (por**m) * ((S**n)* (pore_ec,) + ((por**-m) - 1)*cond_sup)


def mod4linde06(water, bd, pd, clay, pore_ec):
    """Linde et al. 2006 conductivity equation. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below), the effect seems negligible. 'm' is calculated using prop12. Surface conductivity is calculated using 0.7/1000
    """
    por      = 1 - (bd/pd)    
    alpha    = (-0.46 * clay/100) + 0.71
    S        = water/(por)
    m        = 1/alpha
    cond_sup = 0.7/1000
    n = m

    return     (por**m) * ((S**n)* (pore_ec,) + ((por**-m) - 1)*cond_sup)


def mod5linde06(water, bd, pd, clay_cont, pore_ec,n, cond_sup = .7):
    """Linde et al. 2006 conductivity equation. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below), the effect seems negligible. 'm' is calculated using prop13 with a fixed offset==3. Surface conductivity is also fixed as 0.7 m S/m"""
    por      = 1 - bd/pd
    alpha    = (-0.46 * clay_cont/100) + 0.71
    S        = water/(por)
    m        = 1/alpha
    
    return     (por**m) * ((S**n)* (pore_ec,) + ((por**-m) - 1) * cond_sup/1000)


def wundmod(water, bd, clay_cont, org, pore_ec,air_perm, ps, wp, rest_init = 1389, vol_wat_cont_init = 0.07, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):  
    
    """
    Wunderlich et.al 2013 Equation with a modification in the way which 'L' is introduced. Here is calculated using L = (-1/m)+1, where 'm' is calculated using prop12. This modifications does not seem realible cause in Wunderlich model 'L' refers to the shape of pores, not the particles as refers the relationship L = (-1/m)+1.    
    """
    clay = clay_cont/100
    org = org/100
    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    
    por = 1 - bd/pd
    alpha = (-0.46 * clay) + 0.71
    moist = water
    #m = np.log((((1 - por)*(ps/wp)**alpha) + por)) / (np.log(por)*alpha)
    offset = 4
    
    m = (np.log((((1 - por)*(5/80)**alpha) + por)**(1/alpha))-offset/80) / np.log(por)
    L = (-1/m) + 1
    
    dif = moist - vol_wat_cont_init                                       # Diference utilized just for simplicity
    y = rest_init**-1                                                     # Initial permitivity = Epsilon sub 1  
    x = 0                                                                 # Diferentiation from p = 0  
    dx = 0.001                                                            # Diferentiation step
    
    while x<1:                                                            # Diferentiation until p = 1
        dy = dx*((((y)*(dif))/(1-dif+x*dif)) * (((pore_ec,)-(y))/((L/pore_ec) + (1-L)*y)))
        x=x+dx
        y=y+dy
                                                                          # Return electrical conductivity of the soil
    return y


def modwaxmansmits(water, bd, clay_cont, org, CEC, pore_ec,air_perm, ps, wp, n, orgdens = 1.3, B=7.7e-8, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """
    Maxman and Smits 1968 model. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below). 'm' is calculated usinf prop12 with a fixed offset=4 and alpha with the empirical relationship  linkig it with clay content. Exponent for saturation is "n", no simply 2"
    """
    clay = clay_cont/100
    org = org/100
    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    por = 1 - bd/pd
    alpha = (-0.46 * clay) + 0.71
    S = water/por  
    offset = 4
    
    m = (np.log((((1 - por)*(5/80)**alpha) + por)**(1/alpha))-offset/80) / np.log(por)
    Qv = pd*((1-por)/por)*CEC                        # Excess of surface charge per unit volume
    
    return ((por**m)*(S**n)/pore_ec) + ((por**m)*B*Qv/S)


def mod2waxsmidt(water, bd, pdn, CEC, pore_ec,m, n, B0):
    """Waxman & Smidt 1968 modification.
       Exponent for saturation is "n", no simply 2"""
    e = 2.718280828
    B = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por = 1 - bd/pdn
    F = por**-m
    S = water/por                                        # Saturation 
    Qv = pdn*((1-por)/por)*CEC  
            
    return (S**n)/(F*pore_ec)   +  B*Qv*S/F


def mod3waxsmidt(water, bd, clay_cont, org, CEC, pore_ec,air_perm, ps, wp, n, offset, B0, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """
    Maxman and Smits 1968 model. Particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below). 'm' is calculated usinf prop9 where alpha with the empirical relationship  linkig it with clay content. Exponent for saturation is "n", no simply 2"
    """
    clay     = clay_cont/100
    org      = org/100
    somr     = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd       = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    
    por      = 1 - bd/pd
    alpha    = (-0.46 * clay) + 0.71
    S        = water/(por)
    rate     = ps / wp
    m        = 1/alpha
    e = 2.718280828
    B = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por = 1 - bd/pd
    F = por**-m
    S = water/por                                        # Saturation 
    Qv = pd*((1-por)/por)*CEC  
    
    return     (por**m) * ((S**n)* (pore_ec,)) + ((B*Qv)/(F))


def mod4waxsmidt(water, bd, clay_cont, org, CEC, pore_ec,air_perm, wp, m, n, B0, densorg = 1.4, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """Maxman and Smits 1968 model. Particle density as function of organic matter and clay content as is described in Schjonning et al  2012 (def schjonnpd below). Exponent for saturation is "n", no simply 2
    """
    clay     = clay_cont/100
    org      = org/100
    somr     = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd       = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
 
    por      = 1 - bd/pd
    S        = water/(por)
    e        = 2.718280828
    B        = B0*( 1-0.6*np.exp(-1/(0.013*pore_ec)))
    por      = 1 - bd/pd
    F        = por**-m
    S        = water/por                                        # Saturation 
    Qv       = pd*((1-por)/por)*CEC  
    
    return     (por**m) * ((S**n)* (pore_ec,)) + ((B*Qv)/(F))


def modglover10(water, bd, pdn, clay_cont, cond_matrix, pore_ec,ps, wp, n, offset):
    """ Glover (2010b) Derived from the conventional Archie’s law by considering boundary conditions implied by geometric constraints.
    Tabla2, Glover 2015. Three phases. The modification is introduced through 'm' using prop9
    """
    por =    1 - bd/pdn
    alpha = -0.46*clay_cont/100 + 0.71 
    S =     water/por 
    rate =  ps/wp
    m =     1/alpha
    m2 =    np.log(1-water**m)/np.log(1-water)   # This exponent needs also to consider the exponent of the air phase (neglected)
    
    return  (pore_ec,)*water**m + cond_matrix*(1-por)**m2


def LRMod1(water, bd, pdn, clay_cont, cond_matrix, pore_ec):    
    """Lichtenecker and Rother (1931), Korvin (1982) Derived from the theory of functional equations under appropriate
boundary conditions. Formally, the same as Archie’s law if s1=0. Glover 2015 Table 2. 'm' is calculated using prop9"""
    por =   1 - bd/pdn
    alpha = -0.46 * clay_cont/100 + 0.71
    m =     1/ alpha
    
    return  (((pore_ec,)**(1/m))*water + (cond_matrix**(1/m))*(1-por) )**m


def LRMod2(water, bd, pdn, clay_cont, cond_matrix, pore_ec,ps, wp):    
    """Lichtenecker and Rother (1931), Korvin (1982) Derived from the theory of functional equations under appropriate
boundary conditions. Formally, the same as Archie’s law if s1=0. Glover 2015 Table 2. 'm' is calculated using prop11"""
    rate  = ps/wp
    alpha = -0.46 * clay_cont/100 + 0.71
    por   = 1 - bd/pdn     
    
    m     = np.log(((1 - por)*(rate)**alpha) + por) / (np.log(por)*alpha)
    
    return (((pore_ec,)**(1/m))*water + (cond_matrix**(1/m))*(1-por) )**m


def LRMod3(water, bd, pdn, clay_cont, cond_matrix, pore_ec):    
    """Lichtenecker and Rother (1931), Korvin (1982) Derived from the theory of functional equations under appropriate
boundary conditions. Formally, the same as Archie’s law if s1=0. Glover 2015 Table 2.  'm' is calculated using prop13"""
    por =    1 - bd/pdn
    S      = water/por
    
    if [clay_cont >= 5]:
        m = 0.92*clay_cont**0.2
        
    else:
        m = 1.25
        
    return   (((pore_ec,)**(1/m))*water + (cond_matrix**(1/m))*(1-por) )**m


def bardonpiedMod1(water, bd, pdn, clay_cont, cond_clay, pore_ec,ps, n):
    """
    Modification of Bardon and Pied 1969 in Table 4 Glover 2015
    'm' is calculated as in prop11 and the exponent of S is 'n' not just 2"""
    
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    alpha = -0.46 * clay_cont/100 + 0.71
    rate  = ps/80
    offset= 4
    
    m     = (np.log((((1 - por)*(rate)**alpha) + por)**(1/alpha) - (offset/80))) / np.log(por)
    S     = water/por
    F     = por**-m
    
    return ((S**2)/(pore_ec*F)) + clay_cont*cond_clay*S


def schlumberger1(water, bd, pdn, clay_cont, cond_clay, pore_ec):
    """ 
    Modification 1 of schlumberger eq in Table 4 Glover 2015 """
    
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    
    if [clay_cont >= 5]:
        m = 0.92*clay_cont**0.2
    else:
        m = 1.25
        
    F = por**-m
    n = m
    
    return (S**2)/(F*(1-clay_cont)*pore_ec) + clay_cont*cond_clay*S


def schlumberger2(water, bd, pdn, clay_cont, cond_clay, pore_ec,ps, wp, offset):
    """ 
    Modification 2 of schlumberger eq in Table 4 Glover 2015 """
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    rate = ps / wp
    
    alpha = -0.46 * clay_cont/100 + 0.71
    m = (np.log((((1 - por)*(rate)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    F = por**-m
    
    return (S**2)/(F*(1-clay_cont)*pore_ec) + clay_cont*cond_clay*S


def schlumberger3(water, bd, pdn, clay_cont, cond_clay, pore_ec):
    """ 
    Modification 3 of schlumberger eq in Table 4 Glover 2015 """
    clay_cont = clay_cont/100
    por = 1 - bd/pdn
    S = water/por
    alpha = -0.46 * clay_cont/100 + 0.71
    m = 1/alpha
    F = por**-m
    
    return (S**2)/(F*(1-clay_cont)*pore_ec) + clay_cont*cond_clay*S


 ###############################################    REAL SURFACE ELECTRICAL CONDUCTIVITY  ####################################

    
def doussan_ruy_scond(sand, silt, clay):
    " Surface conductivity based on WS model obtained for clay contents higher than 6 %. Based on Rhoades model otherwise. As defined in Doussan & Ruy 2009. R2=0.968"
    if clay > 6:
        return (0.654*clay/(sand+silt)) + 0.018
    else:
        return 0.023*clay - 0.0209
    
    
    
 ###############################################    IMAGINARY ELECTRICAL CONDUCTIVITY  ####################################


def logsdon_eq9i(fc, rperm, iperm, np):
    """
    Logsdon et al 2010 eq 8 for imaginary conductivity
    """
    pi = 3.141592653589793
    epsilon_0 = 8.8541878028e-12
    arc = np.arctan(2.0*iperm/rperm)
    return arc*iperm*pi*fc*epsilon_0


def slater(por, ssa):
    """
    Empirical power law for imaginary conductivity as shown in Glover (2015), taken from Slater (2007)
    """
    spor = ssa/por
    return spor**0.8

    
 ###############################################    Real PERMITIVITY EQUATIONS   ####################################


def MG_2phase(ps, air_perm, bd, pdn):
    """
    The Maxwell-Garnett [1904] mixing model based on the Lord Rayleigh [1892] formula is the most commonly used model for describing a twophase mixture (Robinson & Friedman 2003)"""
    por = 1 - bd/pdn
    f =   1 - por
    
    return air_perm + 3*f*air_perm*((ps - air_perm)/(ps + 2*air_perm - f*(ps - air_perm)))


def bosch_silt(water):
    """Comparison of Capacitance-Based Soil Water Probes in Coastal Plain Soils. David D. Bosch. 2004. bulk density 1.54. sand 87 silt 7 clay 6. Bottom leaching drying experiment. Real permittivity vs moisture content using HydraProbe"""
    p = [ 1.632e-5, -9.751e-4, 3.251e-2, -0.0863 - water ]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    return (perm_rel[0].real)


def bosch_sand(water):
    """Comparison of Capacitance-Based Soil Water Probes in Coastal Plain Soils. David D. Bosch. 2004. bulk density 1.57. sand 90 silt 7 clay 3. Bottom leaching drying experiment. Real permittivity vs moisture content using HydraProbe"""
    p = [ 7.587e-6, -9.331e-4, 3.861e-2, -0.1304 - water ]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    return (perm_rel[0].real)
    
    
def bosch_clay(water):
    """Comparison of Capacitance-Based Soil Water Probes in Coastal Plain Soils. David D. Bosch. 2004. bulk density 1.51. sand 60 silt 9 clay 31. Bottom leaching drying experiment. Real permittivity vs moisture content using HydraProbe"""
    p = [ 3.350e-5, -2.519e-3, 6.625e-2, -0.2093 - water ]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    return (perm_rel[0].real)
    
    
def ojo_MB1(perm):
    """ Laboratory calibrated equation for sandy soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (0.1127* perm**0.5)-0.2025


def ojo_MB4(perm):
    """ Laboratory calibrated equation for sandy soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return perm*(1.085e-2) - (1.2122e-2)*perm**(0.5) + 3.082e-2


def ojo_MB6(perm):
    """ Laboratory calibrated equation for clay soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return  (6.090e-3)*perm + (2.0246e-2)*perm**0.5 + 9.746e-3


def ojo_MB7(perm):
    """ Laboratory calibrated equation for sandy soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (0.1084* perm**0.5)-0.1949


def ojo_MB8(perm):
    """ Laboratory calibrated equation for clay soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return  8.066e-3*perm + 8.76e-4*perm**0.5 + 3.9367e-2


def ojo_MB9(perm):
    """ Laboratory calibrated equation for sandy soil using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (0.1131* perm**0.5)-0.2116


def ojo_all(perm):
    """ Laboratory calibrated equation for all the soils using Hydraprobe. Calibration and Evaluation of a Frequency Domain
Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return 0.0870*perm**0.5 - 0.1425


def ojo_15cec(rperm):
    """ FIELD calibrated equation for soils with less than 15 meq/100gr CEC using Hydraprobe (Table 5, R2=0.94). Calibration and Evaluation of a Frequency Domain Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (0.1084*rperm**(0.5))-0.1633
    
    
def ojo_15_30cec(rperm):
    """ FIELD calibrated equation for soils with 15 to 30 meq/100gr CEC using Hydraprobe (Table 5, R2=0.87). Calibration and Evaluation of a Frequency Domain Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (0.0786*rperm**(0.5)) - 0.0714


def ojo_30cec(rperm):
    """ FIELD calibrated equation for soils with more than 30 meq/100gr CEC using Hydraprobe (Table 5, R2=0.77). Calibration and Evaluation of a Frequency Domain Reflectometry Sensor for Real-Time Soil Moisture Monitoring E. RoTimi Ojo"""
    return (3.32*10**-3)*rperm**(1.5) - (6.784*10**-2)*rperm + (5.047*10**-1)*rperm**(0.5) - 8.85*10**-1


def hydraprobe(water):
    """
        Hydraprobe default equation for water content (See Hydraprobe manual, equation A2 apendix C).

        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    A = 0.109
    B = -0.179
    
    return (((water - B)/A)**2)*epsilon_0


def Hydra_moist(rperm):
    """
    Hydraprobe default equation for VMC (eauqtion A2, apendix C)
    """
    A = 0.109
    B = -0.179
    
    return np.sqrt(rperm)*A + B


def topp(water):
    """
        Topp et al. (1980).
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        Returns
        -------
        Aparent bulk permittivity: float        
    """
    p = [4.3e-6, -5.5e-4, 2.92e-2, -5.3e-2 - water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0


def logsdonperm(water):
    """
        Logsdon 2010. 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        Returns
        -------
        Real bulk permittivity: float            
    """
    
    p = [0.00000514, -0.00047, 0.022, 0-water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0
    
    
def nadler(water):
    """
        Nadler et al. (1991). 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        Returns
        -------
        Aparent bulk permittivity: float            
    """
    
    p = [15e-6, -12.3e-4, 3.67e-2, -7.25e-2 - water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0


def roth1992(water):
    """
    Roth et al. (1992) 
    """
    p = [36.1e-6, -19.5e-4, 4.48e-2, -7.28e-2 - water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0


def jacandschjA(water):
    """
        Jacobsen and Schjonning (1993) (Equation 1)
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]  
            
        Returns
        -------
        Aparent bulk permittivity: float            
    """
    
    p = [18e-6, -11.6e-4, 3.47e-2, -7.01e-2 - water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0


def jacandschjB(water, bd, cc, org):
    """
    Jacobsen and Schjonning (1993) (Equation 2)
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        cc: float
            volumetric clay content (%)
            
        org: float
            organic matter content (%)
            
        Returns
        -------
        Aparent bulk permittivity: float
    """
    
    p = [17.1e-6, -11.4e-4, 3.45e-2, 
         -3.41e-2 - water -3.7e-2 * bd + 7.36e-4 * cc + 47.7e-4 * org]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    
    return (perm_rel[0].real)*epsilon_0


def hallikainen_1_4(water, cc, sand):
    """ 
        Empirical model for permittivity at 1.4 Ghz. Hallikainen et al., 1985
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        cc: float
            volumetric clay content (%)
            
        sand: float
            volumetric sand content (%)
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    a0 = 2.862 
    a1 = -0.012 
    a2 = 0.001 
    b0 = 3.803 
    b1 = 0.462 
    b2 = -0.341 
    c0 = 119.006 
    c1 = -0.500 
    c2 = 0.633 
    
    return ((a0 + a1*sand + a2*cc) + (b0 + b1*sand + b2*cc)*water + (c0 + c1*sand + c2*cc)*water**2)*epsilon_0


def hallikainen_4(water, cc, sand):
    """ 
        Empirical model for permittivity at 4 Ghz. Hallikainen et al., 1985
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        cc: float
            volumetric clay content (%)
            
        sand: float
            volumetric sand content (%)
            
        Returns
        -------
        Real bulk permittivity: float            
    """
    
    a0 = 2.927
    a1 = -0.012 
    a2 = -0.001 
    b0 = 5.505 
    b1 = 0.371 
    b2 = 0.062
    c0 = 114.826
    c1 = -0.389
    c2 = -0.547
    
    return ((a0 + a1*sand + a2*cc) + (b0 + b1*sand + b2*cc)*water + (c0 + c1*sand + c2*cc)*water**2)*epsilon_0

    
def hallikainen_18(water, cc, sand):
    """ 
        Empirical model for permittivity at 18 Ghz. Hallikainen et al., 1985
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        cc: float
            volumetric cc content [-]
            
        sand: float
            volumetric sand content [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    a0 = 1.912
    a1 = 0.007
    a2 = 0.021
    b0 = 29.123 
    b1 = -0.190 
    b2 = -0.545
    c0 = 6.960
    c1 = 0.822
    c2 = 1.195
    
    return ((a0 + a1*sand + a2*cc) + (b0 + b1*sand + b2*cc)*water + (c0 + c1*sand + c2*cc)*water**2)*epsilon_0


def raizfunc(water):

    return(((water + 0.1788)/0.1138)**2)*epsilon_0 


def steelman(water):
    """
        Steelman* and Anthony L. Endres (2011) 
              
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
            
        Returns
        -------
        Aparent bulk permittivity: float            
    """
    
    p = [2.97e-5, -2.03e-3, 5.65e-2, -0.157 - water]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    perm_rel = roots[roots > 0]
    return (perm_rel[0].real)*epsilon_0


def malicki(water, bd):
    """
        Malicki et al. 1996
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        bd: float
            bulk density [g/cm3]
            
        Returns
        -------
        Aparent bulk permittivity: float
    """
    
    return((water*(7.17 + 1.18*bd) + 0.809 + 0.168*bd + 0.159*bd**2)**2)*epsilon_0
    
    
def crim(water, bd, pdn, ap, sp, wp, alpha = 0.5):
    """
        Birchak et.al 1974
        
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
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn        # Porosity
    
    return (( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha))*epsilon_0


def roth1990_bw(water, bwc, bd, pdn, air_perm, ps, wp, bwp, alpha):
    """
    # Roth et.al (1990) Equation
    """
    por = 1 - bd/pdn
    
    return (( water*wp**alpha + bwc*bwp**alpha + (1-por)*ps**alpha + (por-water)*air_perm**(alpha))**(1/alpha))*epsilon_0


def roth(water, bd, pdn, ap, sp, wp, alpha): 
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
            alpha exponent [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn    # Porosity
    
    return (( water*wp**alpha + (1-por)*sp**alpha + (por-water)*ap**(alpha))**(1/alpha))*epsilon_0


def linde_p(water, bd, pdn, ap, sp, wp, m, n):
    """
        Linde et al., 2006
        
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
            
        m: float
            cementation exponent/factor [-]
            
        n: float
            saturation exponent [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn  # Porosity
    S = water / por     # Saturation
    
    return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*sp)+(1-S**n)*ap)*epsilon_0


def sen81eq23(bd, pdn, ps, wp):
    """ Sen 1981 equation 23 for saturated conditions and spherical grains, in DC limit """
    por = 1 - bd/pdn
    
    return (1.5*ps + (por**1.5)*(wp - 1.5*ps))*epsilon_0


def pride94(bd, pdn, ps, wp, m):
    """
    Pride 1994 Equation for permittivity (eq. 253)
    """
    por = 1 - bd/pdn
    
    return ((por**m) * (wp - ps) + ps)*epsilon_0


def cornelis(water, bd, pdn, ps):
    """
    Equation shown in Corneli's course
    """
    return ((1 + (((ps**0.5) - 1) *  bd)/pdn + 8*water)**2)*epsilon_0


def dobson(water, bw, bd, pdn, ap, sp, wp, bwp):
    """
        Dobson et al., 1985  

        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        bw: float
            volumetric bound water content [-]
            
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
            
        bwp: float
            bound water permittivity [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn
    
    num = 3*sp+2*water*(wp-sp) + 2*bw*(bwp-sp)+2*(por-water)*(ap-sp)
    den = 3+water*((sp/wp)-1) + bw*((sp/bwp)-1)+(por-water)*((sp/ap)-1)
    
    return (num/den)*epsilon_0



def sen(water, bd, pdn, ap, sp, wp, L):        
    """
        Sen et al., 1981
        
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
            
        L: float
            depolarization factor [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn                                       # Porosity 
    cl = water*(L*sp + (1-L)*wp)                             # Calculation just for simplicity 
    wcg = wp*(((1-por)*sp+cl) / ((1-por)*wp+cl))           # water coated grains
    df = (por*-1) + water                                    # Diference utilized just for simplicity
                                                           # Initializing diferenciation parameters
    y = ap                                                 # Initial permitivity = epsilon sub a  
    x = 0.001                                              # Diferentiation from p = 0  
    dx = 0.05                                              # Diferentiation step
                                                           
    while x<1:                                             # Diferentiation until p = 1
        dy = ((y*(1+df))/(-df+x*(1+df))) * ((wcg-y)/(L*wcg+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y*epsilon_0 
 
    
def feng_sen(water, bd, pdn, ap, sp, wp, L):
    """
        Feng & Sen 1985
        
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
            
        L: float
            depolarization factor [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn                                       # Porosity
                                                           # Initializing diferenciation parameters
    y = wp                                                 # Initial permitivity = epsilon sub a  
    x = 0                                                  # Diferentiation from p = 0  
    dx = 0.05                                              # Diferentiation step
    
    while x<1:                                             # Diferentiation until p = 1
        dy = (y/(water+x*(1-water))) * ((((1-por)*((sp-y))/(L*sp+(1-L)*y))) + ((por-water)*(ap-y))/(L*ap+(1-L)*y)) 
        x = x + dx
        y = y + dy*dx
        
    return y*epsilon_0 



def endres_redman(water, bd, pdn, ap, sp, wp, L):   
    """
        Endres & Redman 1996
        
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
            
        L: float
            depolarization factor [-]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn                                              # Porosity
    S = water/por                                                   # Saturation
                                                                  # Initializing diferenciation parameters
    y = wp                                                        # Initial permitivity = epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.05                                                     # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = ((dx*y*(1-S))/(S+x*(1-S))) * ((ap-y)/(L*ap+(1-L)*y))  
        x = x + dx
        y = y + dy
                                                                  # Now y is equal to permitivity of pore(s)
    p = 0
    dp = 0.05
    z = y
    
    while p<1:    
        dz = (dp*z*(1-por))/(por+p*(1-por)) * ((sp-z)/(L*sp+(1-L)*z))
        p = p + dp
        z = z + dz
        
    return z*epsilon_0




def peplinski(water, bd, pdn, cc, sand, sp, wp, ewinf, tau = 0.65):
    """
        Peplinski et al., 1995 (Equations 1 to 6) 
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]
         
        cc: float
            volumetric clay content (%)
           
        sand: float
            volumetric sand content (%)
            
        sp: float
            solid permittivity phase [-]
            
        wp: float
            water permittivity phase [-]
            
        ewinf: float
            permittivity at infinite frequency [-]
            
        Returns
        -------
        Real bulk permittivity: float   
    """
    
    beta1 = 1.2748 - 0.519*sand/100 - 0.152*cc/100
    
    return (1 + (bd/pdn)*sp**tau + (water**beta1)*(wp**tau) - water)**(1/tau)*epsilon_0


def stratton41(bd, pdn, L, ps, wp):
    """
    Eq A-1 Sen et al. 1981
    """
    por = 1 - (bd/pdn)
    bulkperm = wp*((wp+(ps-wp)*(1-por+por*L))/(wp + por*L*(ps-wp)))
    return bulkperm*epsilon_0


def strattonL0(water, bd, pdn, ps, wp):
    """Sen 1981 Eq A-3 derived using L=0 in Stratton 1941 equation. It is valid for complex permittivity. 
    Porosity along with wp was replaced for moisture content."""
    por = 1- bd/pdn
    return ((1 - por)*ps + water*wp)*epsilon_0


def strattonL1(water, bd, pdn, ps, wp):
    """Sen 1981 Eq A-4 derived using L=1 in Stratton 1941 equation. It is valid for complex permittivity.
    Porosity along with wp was replaced for moisture content."""
    por = 1- bd/pdn
    return (((1 - por)/ps + water/wp)**-1)*epsilon_0


def logsdon(fq, rperm, iperm):
    """
        Equation 8 of Logsdon 2010
    
        Parameters
        ----------
        fq: float
            working frequency (Hz)
            
        rperm: float
            real permittivity
        
        iperm: float
            imaginary permittivty
            
        Returns
        -------
        Real bulk direct current electrical conductivity: float   
    """
    
    arc = np.arctan(2.0*iperm/rperm)
    return arc*rperm*pi*fq*epsilon_0


def debye(freq, time, eps_0, eps_inf):
    """Debye model for frequency dependence of polarization process as defined in Glover 2015 pag 27. The perm utilized can be refered to a mean values for different moisture contents."""
    """ The value for ε∞ = 4.9, while εfws and τfw are calculated as a function of pore water salinity and temperature following Stogryn (1971)."""
    return eps_inf + (eps_0-eps_inf)/(1+((2*np.pi*freq)**2)*time*np.exp2)
    

def cole_cole(freq, eps_0, relax_time=2e-7, alpha_cole=0.7, eps_inf=4.9):
    """A modification of the Debye model by Cole and Cole (1941) takes an ensemble of relaxation times into account(Glover 2015)"""
    """ The value for ε∞ = 4.9, while εfws and τfw are calculated as a function of pore water salinity and temperature following Stogryn (1971)."""
    """Alpha_cole and relax_time is taken from: 'Estimation of frequency domain soilparameters of horizontally multilayered earthby using Cole–Cole model based on theparallel genetic algorithm'  """
    
    return eps_inf + (eps_0-eps_inf)/ (1+(2*np.pi*freq*relax_time)**(1-alpha_cole))


def MWBH(    bd, pd, rsw, ps, freq, wp):
    """"Permittivity WMBH model for frequency between 100Mhx and 1GHz as described (eq 5) and tested in Chen, Y., and D. Or (2006). Here the samples are sands and sandy loam almost saturated, while debye is not."""
    epsilon_0 =      8.8541878028e-12                                   # Vacum permittivity
    pi =       3.141592653589793
    por =      1 - bd/pd
    pors =     1 - por
    deltae =   (9**pors*(1-pors)*ps**2) / (((2+pors)**2)*(2*wp+ps+pors*(wp-ps)))
    fr =       (2+pors)/(rsw*2*np.pi*epsilon_0*(2*wp+ps+pors*(wp-ps)))  # Relaxation frequency
    t =        1/(2*np.pi*fr)
    sigma =    2*(1-pors)/(rsw*(2+pors))
    w    =     2*np.pi*freq
    perm =     (wp*(2*(1-pors)*wp + (1+2*pors)*ps) / (wp*(2+pors)+ps*(1-pors))) + deltae/(1+ (w*t)**2) + sigma/(w*epsilon_0)
    
    return perm
    
    
################################ REAL Water PERMITTIVITY VS TEMPERATURE #########################################


def jones05(T):
    """ Jones 2005 equation as defined in HydraProbe manual for water permittivity"""
    #T = T + 273.15
    perm = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    return perm


def handbook(T):
    """ Weast, R.C. 1986. CRC Handbook of Chemistry and Physics. CRC
soils. Applied Geophysics 49:73–88. Press, Boca Raton, FL."""
    return 78.5411- (4.579e-3)*(T-25)+ (1.19e-5)*(T-25)**2 - (2.8e-8)*(T-25)**3


def teros12pt(T):
    """ Water permittivity vs T as defined in eq 3 in Teros12 manual """
    return 80.3 - 0.37*(T-20)


##################################    REAL BULK PERMITTIVITY MODELS FOR Temperature #############################


# Bichark et al. (1974) shows a combination of volumetric mixing model where the temperature is introduced trough the water phase permittivity. This is the base of all the models described in this section. See for instance Seyfried & Murdock 2004 page 401.
def tCRIM(water, bd, pdn, ps, T, air_perm = 1, alpha = 0.5):
    """
    #CRIM Equation (Birchak et.al 1974)
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    
    return ((water * ((wp**alpha) - (air_perm**alpha)) + ((1-por)*(ps**alpha)) + por*(air_perm**alpha))**(1/alpha))


def troth1990(water, bd, pdn, ps, T, air_perm = 1): 
    """
    # Roth et.al (1990) Equation
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    alpha = 0.5    
    
    return ((water * ((wp**alpha) - (air_perm**alpha)) + ((1-por)*(ps**alpha)) + por*(air_perm**alpha))**(1/alpha))


def tsen1981(water, bd, pdn, L, ps, T, air_perm = 1):        
    """
    # Sen et.al 1981 Equation
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    cl = water*(L*ps + (1-L)*wp)                            # Calculation just for simplicity 
    wcg = wp*(((1-por)*ps+cl) / ((1-por)*wp+cl))   # wcg = wat coated grains
    df = (por*-1) + water                                          
    y = air_perm                                                                 
    x = 0.001                                                                        # Diferentiation from p = 0  
    dx = 0.01                                                                        # Diferentiation step
    
    while x<1:                                                                       # Diferentiation until p = 1
        dy = ((y*(1+df))/(-df+x*(1+df))) * ((wcg-y)/(L*wcg+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y
   
    
def tfengsen1985(water, bd, pdn, L, ps, T, air_perm = 1):
    """
    # Feng & Sen 1985 Equation
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    water = water                                          # Abreviation just for simplicity 
    pds = ps                                          # Abreviation just for simplicity
                                                                  # Initializing diferenciation parameters
    y = wp                                                # Initial permitivity = Epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.001                                                     # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = (y/(water+x*(1-water))) * ((((1-por)*((pds-y))/(L*pds+(1-L)*y))) + ((por-water)*(air_perm-y))/(L*air_perm+(1-L)*y)) 
        x = x + dx
        y = y + dy*dx
        
    return y


def tendresredman1996(water, bd, pdn, L, ps, T, air_perm = 1):   
    """
    # Endres & Redman 1996 Equation
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    water = water                                          # Abreviation just for simplicity 
    pds = ps                                          # Abreviation just for simplicity
    S = water/por                                              # Saturation
                                                                  # Initializing diferenciation parameters
    y = wp                                                # Initial permitivity = Epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.001                                                     # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = ((dx*y*(1-S))/(S+x*(1-S))) * ((air_perm-y)/(L*air_perm+(1-L)*y))  
        x = x + dx
        y = y + dy
                                                                  # Now y is equal to permitivity of pore(s)
    p = 0
    dp = 0.001
    z = y
    
    while p<1:    
        dz = (dp*z*(1-por))/(por+p*(1-por)) * ((pds-z)/(L*pds+(1-L)*z))
        p = p + dp
        z = z + dz
        
    return z


def twunderlich2013perm(water, L, T, vol_wat_cont_init = 0.04, perm_init = 6): 
    """
    #Wunderlich et.al 2013 Equation for permitivity
    """
    # Taking into account the fact that clay content has an appreciable influence in wp. Pag 10/14 Wunderlich et al.2013
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)

    diff = water - vol_wat_cont_init                      # Diference utilized just for simplicity
    y = perm_init                                                  # Initial permitivity = Epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.001                                                     # Diferentiation step
                                                                   # Diferentiation until p = 1
    while x<1:                                                    
        dy = ((y*diff)/(1-diff+x*diff)) * ((wp-y)/(L*wp+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y


def tpride94(water, bd, pdn, L, ps, n, T):
    """
    Pride 1994 Equation for permittivity (eq. 253)
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    m = 1/(1-L)
    
    return ((por**m) * (wp - ps) + ps)


def tlindeperm(water, bd, pdn, L, ps, n, T, air_perm = 1):
    """
    Linde permittivity 2006
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    por = 1 - bd/pdn
    m = 1/(1-L)
    S = water / por
    #The folowing equation takes intp account the influence of air in bulk permittivity, the changes are minimum
    #return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*ps + (1 - S**n)*air_perm))*epsilon_0
    return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*ps))


def tpeplinski(water, bd, pdn, clay_cont, sand_cont, ps, frec, T, ewinf, alpha = 0.65, tau = 0.58e-10):
    """
    Peplinski 1995 model for real component of permittivity eqs 1 - 6 
    """
    wp = 80*(1-(4.579e-3) *(T-25)+ (1.19e-5)* ((T-25)**2) - (2.8e-8)*(T-25)**3)
    robulkrosolid = bd/pdn
    beta1 = 1.2748 - 0.519*sand_cont/100 - 0.152*clay_cont/100
    efw = ewinf + (wp - ewinf)/(1 + (tau*frec)**2)
    
    bulkperm = (1 + robulkrosolid*ps**alpha + (water**beta1)*(efw**alpha) - water)**(1/alpha)
    return bulkperm


##########################################   SOLID PHASE PERMITTIVITY   ###########################################


def crim_es(bulkperm, bd, pdn, ap):
    """
        Kameyama & Miyamoto, 2008 
        
        Parameters
        ----------
        bulkperm: float
            bulk real permittivity [-]
        
        bd: float
            bulk density [g/cm3]
        
        pdn: float
            particle density [g/cm3]
            
        ap: float
            air permittivity phase [-]

        Returns
        -------
        Solid permittivity phase: float   
    """
    
    por = 1 - bd/pdn                                           # Porosity
    return ((bulkperm**0.5 - por*ap**0.5)/(1-por))**2


def linde_es(bulkperm, bd, pdn, m):
    """
    Linde et al. 2006 model for solid phase permittivity"""
    por = 1 - bd/pdn
    F = por**-m
    return (F/(F-1))*bulkperm


def olhoeft_b(bulkperm, bd, pdn):
    """
    'Olhoeft [1981] used the expression (1a) based on the Lichtenecker [1926] equation (1b), which simply averages the logarithms of the permittivities' (Robinson & Friedman 2003)"""
    por = 1 - bd/pdn
    return bulkperm**(1/(1-por))
    
    
def nelson(bulkperm, bd, pdn):
    """
    'Nelson et al. [1989] favored the Looyenga [1965] mixing formula' (Robinson & Friedman 2003)"""
    por = 1 - bd/pdn
    f = 1 - por
    return ((bulkperm**0.3 + f - 1)/f)**3


def bruggeman(bulkperm, air_perm, bd, pdn):
    """
    Bruggeman formula, the symmetric effective medium approximation [Bruggeman, 1935] (Robinson & Friedman 2003)"""
    por = 1 - bd/pdn
    F =   (-por*(air_perm**0.5)) + 1
    return (1 - F - bulkperm**(2/3)) / (1 - F - bulkperm**(-1/3))


def Dobson(part_dens):
    """" Ca is 1.0 and es is determined by an empirical fitting of the data presented in Part I for soils having extremely low moisture contents. The resultant expression es = (1.01 + 0.44 p,)2 - 0.062 (22) (at frequencies between 1.4 and 18 GHz.) yields es 4.7 at the specific densities given in Table 1. Equation(22) is nearly identical to that obtained by Shutko [21] and Krotikov [22] for other soils.(Microwave Dielectric Behavior of Wet Soil-Part II: Dielectric Mixing Models, Dobson 1985)"""
    return (1.01+0.44*part_dens)**2-0.062


def HBS_es(bd, pdn, bulk_perm, air_perm):
    "Bruggeman 1935; Hanai 1968; Sen et al. 1981 model for a mix of spherical grains and air (Guillemoteau etal 2011)"
    por = 1 - bd/pdn
    const = por*(air_perm/bulk_perm)**(-1/3)  # Just for simplify writing
    
    return (-const*air_perm + bulk_perm)/(1 - const)
    
    
##########################################   Imaginary PERMITTIVITY     ###########################################
    
    
def ipeplinski(water, bd, pdn, clay_cont, sand_cont, frec, wp, alpha = 0.65 , ewinf = 4.9, tau = 0.58e-10):
    """
    Peplinski 1995 model for imaginary component of permittivity eqs 1 - 6 
    """
    beta2 = 1.33797 - 0.603*sand_cont/100 - 0.166*clay_cont/100
    cond = 0.0467 + 0.2204*bd - 0.4111*sand_cont/100 + 0.6614*clay_cont/100
    
    efw2 = ((tau*frec*(wp - ewinf))/(1+(tau*frec)**2)) + cond*(pdn - bd)/(2*pi*epsilon_0*frec*pdn*water)
    
    iperm = ((water**beta2)*(efw2**alpha))**(1/alpha)
    
    return iperm*epsilon_0 


def ipeplogsdon(water, bd, pdn, clay_cont, sand_cont, frec, boundwat, wp, alpha = 0.65, ewinf = 4.9, tau = 0.58e-10):
    """
    Logsdon 2010 eq for relaxation component (measured with 50*10**6 Hz frec) and Peplinski 1995 for conductivity component of imaginary permittivity    
    """
    beta2 = 1.33797 - 0.603*sand_cont/100 - 0.166*clay_cont/100
    cond = 0.0467 + 0.2204*bd - 0.4111*sand_cont/100 + 0.6614*clay_cont/100
    perm_eff = -6.55 + 22.3*(water-boundwat) + 161*boundwat
    
    efw2 = perm_eff + (cond*(pdn - bd)*(water**beta2/alpha))/(2*pi*epsilon_0*frec*pdn*water)
    return efw2*epsilon_0 


def sab_ech2(water):
    """ Sabouroux & Ba, Progress In Electromagnetics Research B, Vol. 29, 191{207, 2011. ECH2 sample (pure sand). Pag 7, eq 9."""
    return 0.06 + 1.35*water - 0.53*water**2 + 18.58*water**3


def sab_ech3(water):
    """ Sabouroux & Ba, Progress In Electromagnetics Research B, Vol. 29, 191{207, 2011. ECH3 sample (sand with 10% clay). Pag 7, eq 10."""
    return 0.03 + 4.1*water - 11.3*water**2 + 65.35*water**3
    
    
def sab_ech4(water):
    """ Sabouroux & Ba, Progress In Electromagnetics Research B, Vol. 29, 191{207, 2011. ECH4 sample (sand with 20% clay). Pag 7, eq 11."""
    return 0.14 + 4.63*water - 28.04*water**2 + 151.67*water**3
    
    
def idebye(frec, time, eps_0, eps_inf):
    """Debye model for frequency dependence of polarization process as defined in Glover 2015 pag 27."""
    return (((eps_0-eps_inf)*2*np.pi*frec*time)/(1+((2*np.pi*frec)**2)*time*np.exp2))*epsilon_0


def ieq67(frec, cond, iperm):
    """ Imaginary part of total permittivity as defined in eq 68 Glover etal 2015"""
    return (iperm + cond/(2*np.pi*frec))*epsilon_0


################################################### MODIFIED PERMITTIVITY MODELS ###################################


def modpride94(water, bd, pdn, clay_cont, ps, wp, offset):
    """
    MODIFICATION OF Pride 1994 Equation for permittivity
    """
    por = 1 - bd/pdn
    alpha = (-0.46 * clay_cont/100) + 0.71
    offset = 5
    #m = np.log(((1 - por)*(ps/wp)**alpha) + por) / (np.log(por)*alpha)
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    #m = 1/alpha
    S = water / por
    
    return ((por**m) * (wp - ps) + ps)*epsilon_0


def modroth1990(water, bd, pd, clay_cont, air_perm, ps, wp): 
    """
    MODIFICATION OF Roth et.al (1990) Equation using Wunderlich et al. empirical relationship among roth alpha exponent and clay content. Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below)  
    """
    por = 1 - bd/pd
    alpha = -0.46*clay_cont/100 + 0.71  
    
    return (( water*wp**alpha + (1-por)*ps**alpha + (por-water)*air_perm**(alpha))**(1/alpha))*epsilon_0


def roth_moist(bulk_perm, bd, pdn, air_perm, ps, wp): 
    """
    # Roth et.al (1990) Equation modified with IOP data obtained at 50MHz
    """
    por = 1 - bd/pdn    
    alpha = 0.5
    return ((bulk_perm**alpha)-((1-por)*ps**alpha)-(por*air_perm**alpha))/((wp**alpha)-(air_perm**alpha))


def lindep_moist(bulk_perm, bd, pdn, air_perm, ps, wp, m):
    """
    Linde permittivity 2006 modified for moisture prediction
    """

    por = 1 - bd/pdn
    
    return por * np.exp(1/m*(np.log(   (bulk_perm+(por**m-1)*ps)/(m*wp*por**m)  )))


def wund13permmod(water, bd, clay_cont, org, ps, wp, densorg = 1.3, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209, vol_wat_cont_init = 0.04, perm_init=6): 
    """
    #Wunderlich et.al 2013 Equation for permitivity
    Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below). This function has the issue of calculing L based on 'm' but that L is for particles, not for pores as is defined in Wunderlich model, reason why this is not theoretically viable. 
    """
    clay = clay_cont/100
    org = org/100
    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    
    por = 1 - bd/pd
    alpha = (-0.46 * clay) + 0.71
    offset = 5
    #m = np.log(((1 - por)*(ps/wp)**alpha) + por) / (np.log(por)*alpha)
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    print()
    #m = 1/alpha
    L = (-1/m) + 1
        
    diff = water - vol_wat_cont_init                      # Diference utilized just for simplicity
    y = perm_init                                                  # Initial permitivity = Epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.001                                                     # Diferentiation step
                                                                   # Diferentiation until p = 1
    while x<1:                                                    
        dy = ((y*diff)/(1-diff+x*diff)) * ((wp-y)/(L*wp+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y*epsilon_0



def modlindeperm(water, bd, pd, clay_cont, air_perm, ps, wp, offset):
    """
    Linde permittivity 2006 modification introducing pro3 equation
    Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below)
    """    
    por = 1 - bd/pd
    alpha = -0.46 * clay_cont/100 + 0.71
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    n = m
    S = water / por
    
    return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*ps + (1 - S**n)*air_perm))*epsilon_0


def mod2lindeperm(water, bd, pd, clay_cont, air_perm, ps, wp):
    """
    Linde permittivity 2006 modification introducing pro3 equation
    Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below)
    """    
    por = 1 - bd/pd
    
    if [clay_cont >= 5]:
        m = 0.92*clay_cont**0.2
    else:
        m = 1.25

    n = m
    #m = 1/alpha
    S = water / por
    
    return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*ps + (1 - S**n)*air_perm))*epsilon_0


def modsen1981(water, bd, pd, clay_cont, air_perm, ps, wp, offset):        
    """
    # Sen et.al 1981 Equation modified introducing prop3 link between L and Clya_cont, por, ps and wp
    Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below)
    """
    por = 1 - bd/pd
    
    alpha = -0.46 * clay_cont/100 + 0.71
    #m = np.log(((1 - por)*(ps/wp)**alpha) + por) / (np.log(por)*alpha)
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    #m = 1/alpha
    L = (-1/m) + 1
    #print("Mod Sen81 L = ", L)
    cl = water*(L*ps + (1-L)*wp)                         # Calculation just for simplicity 
    wcg = wp*(((1-por)*ps+cl) / ((1-por)*wp+cl))                  # wcg = wat coated grains
    df = (por*-1) + water                                        # Diference utilized just for simplicity
                                                                          # Initializing diferenciation parameters
    y = air_perm                                                          # Initial permitivity = Epsilon sub a  
    x = 0.001                                                             # Diferentiation from p = 0  
    dx = 0.01                                                             # Diferentiation step
    
    while x<1:                                                            # Diferentiation until p = 1
        dy = ((y*(1+df))/(-df+x*(1+df))) * ((wcg-y)/(L*wcg+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y*epsilon_0 
     
    
def modfengsen1985(water, bd, pd, clay_cont, air_perm, ps, wp, offset, densorg = 1.3, denspart = 2.65, densclay = 2.86, a = 1.127, b = 0.373, c = 2.648, d = 0.209):
    """
    # Feng & Sen 1985 Equation modified introducing prop3 link between L and Clya_cont, por, ps and wp
    Also particle density as function of organic matter and clay content as is described in Schjonning et al 2012 (def schjonnpd below)
    """
    por = 1 - bd/pd
    alpha = -0.46 * clay_cont/100 + 0.71
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    L = (-1/m) + 1
    water = water                                          # Abreviation just for simplicity 
    pds = ps                                              # Abreviation just for simplicity
                                                                  # Initializing diferenciation parameters
    y = wp                                                # Initial permitivity = Epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.001                                                    # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = (y/(water+x*(1-water))) * ((((1-por)*((pds-y))/(L*pds+(1-L)*y))) + ((por-water)*(air_perm-y))/(L*air_perm+(1-L)*y)) 
        x = x + dx
        y = y + dy*dx
        
    return y*epsilon_0  


def linde_mv(water, bd, pdn, ap, sp, wp, CEC):
    """
        Linde et al., 2006 and Mendoza Veirana et al., 2022
        
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
            Cation exchange Capacity [meq/100g]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    m = -0.266*np.log(CEC) + 1.718
    n = m
    por = 1 - bd/pdn
    S = water / por
    
    return ((por**m) * ((S**n)*wp + ((por**-m) - 1)*sp)+(1-S**n)*ap)*epsilon_0    


def sen_mv(water, bd, pdn, ap, sp, wp, CEC):        
    """
        Sen et al., 1981
        
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
            Cation exchange Capacity [meq/100g]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pdn

    L = -0.174 *np.log(CEC) +  0.53
    if CEC> 25:
        L = -0.02
        
    cl = water*(L*sp + (1-L)*wp)                             # Calculation just for simplicity 
    wcg = wp*(((1-por)*sp+cl) / ((1-por)*wp+cl))           # wcg = wat coated grains
    df = (por*-1) + water                                    # Diference utilized just for simplicity
                                                           # Initializing diferenciation parameters
    y = ap                                                 # Initial permitivity = epsilon sub a  
    x = 0.001                                              # Diferentiation from p = 0  
    dx = 0.05       
                                                           # Diferentiation step
    while x<1:                                             # Diferentiation until p = 1
        dy = ((y*(1+df))/(-df+x*(1+df))) * ((wcg-y)/(L*wcg+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y*epsilon_0 


def feng_sen_mv(water, bd, pd, ap, sp, wp, CEC):
    """
        Feng & Sen 1985 and Mendoza Veirana et al., 2022
        
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
            Cation exchange Capacity [meq/100g]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pd

    L = -0.171*np.log(CEC) +  0.378
                                                                  # Initializing diferenciation parameters
    y = wp                                                        # Initial permitivity = epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.05                                                     # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = (y/(water+x*(1-water))) * ((((1-por)*((sp-y))/(L*sp+(1-L)*y))) + ((por-water)*(ap-y))/(L*ap+(1-L)*y)) 
        x = x + dx
        y = y + dy*dx
        
    return y*epsilon_0  


def endres_redman_mv(water, bd, pd, ap, sp, wp, CEC):   
    """
        Endres & Redman 1996
        
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
            Cation exchange Capacity [meq/100g]
            
        Returns
        -------
        Real bulk permittivity: float
    """
    
    por = 1 - bd/pd
    S = water/por          
    L = -0.174 *np.log(CEC) +  0.42
    y = wp                                                        # Initial permitivity = epsilon sub a  
    x = 0                                                         # Diferentiation from p = 0  
    dx = 0.05                                                     # Diferentiation step
    
    while x<1:                                                    # Diferentiation until p = 1
        dy = ((dx*y*(1-S))/(S+x*(1-S))) * ((ap-y)/(L*ap+(1-L)*y))  
        x = x + dx
        y = y + dy
                                                                  # Now y is equal to permitivity of pore(s)
    p = 0
    dp = 0.001
    z = y
    
    while p<1:    
        dz = (dp*z*(1-por))/(por+p*(1-por)) * ((sp-z)/(L*sp+(1-L)*z))
        p = p + dp
        z = z + dz
        
    return z*epsilon_0


def wunderlich_mv(water, perm_init, water_init, wp, CEC): 
    """
        Wunderlich et.al 2013 and Mendoza Veirana et al., 2022
        
        Parameters
        ----------
        water: float
            volumetric moisture content [-]
        
        perm_init: float
            lowest real permittivity [-]
            
        wat_init: float
            lowest volumetric water content [-]
            
        wp: float
            water permittivity phase [-]
            
        CEC: float
            Cation exchange Capacity [meq/100g]
   
        Returns
        -------
        Real bulk permittivity: float   
    """ 
    
    L =  -0.0457  *np.log(CEC) + 0.1181
    diff = water - water_init                                        # Diference utilized just for simplicity
    y = perm_init                                                  # Initial permitivity = epsilon sub 1  
    x = 0.001                                                      # Diferentiation from p = 0  
    dx = 0.05                                                      # Diferentiation step
                                                                   # Diferentiation step until p = 1
    while x<1:                                                    
        dy = ((y*diff)/(1-diff+x*diff)) * ((wp-y)/(L*wp+(1-L)*y))
        x=x+dx
        y=y+dy*dx
        
    return y*epsilon_0


################################### Susc. magnetic ##################################################


def IOPsusc1(Pb, Cu):
    """
    For F1mass"""
    return 0.5949962*Pb + 0.28113054*Cu - 2.1920528397817787


def IOPsusc2(Pb, Cu):
    """
    For F3mass """
    return 0.56453374*Pb + 0.26369148*Cu - 2.095220916590989