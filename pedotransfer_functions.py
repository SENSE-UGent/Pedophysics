"""
    ...

    :AUTHOR: Gaston Mendoza Veirana
    :CONTACT: gaston.mendozaveirana@ugent.be

    :REQUIRES: numpy, scipy
"""

# Import
import numpy as np
from scipy.constants import epsilon_0, pi # TODO replace by numpy constants
e = 2.718280828


##################### Bulk density prediction pedotransfer functions ####################################


def world_wosis_bd(silt, clay, OC, depth):
    """
        Mendoza Veirana et al., 2023

        Parameters
        ----------

        clay: float
            Soil volumetric clay content [m**3/m**3]

        silt: float
            Soil volumetric silt content [m**3/m**3]

        OC: float
            organic carbon [g/kg]

        depth: float
            Soil depth [m]
        
        Returns
        -------
        bulk density [g/cm3]       
    """    

    return 1.7266445721420223 + 6.44286209e-03*clay -9.90540563e-03*OC -3.46590526e-03*silt + 5.20359306e-04*depth**2 -7.26756929e-05*clay**2 + 1.97904579e-05*OC**2 + 3.06762490e-05*clay*OC -2.21456073e-04*depth*clay -1.09073048e-03*depth*OC + 6.08577414e-04*depth*silt


def de_vos(LOI):
    """
        De Vos et al. 2005 
        n = 1614
        R2 = 0.57
        Parameters
        ----------
        LOI: float
            loss on ignition [g/kg]
            
        Returns
        -------
        Bulk density: float    
    """

    return 1.775 - 0.173*LOI**0.5


def abdelbaki(OC):
    """
        Abdelbaki 2018
        n = 45195
        RMSE = 0.14, R2 = 0.68

        Parameters
        ----------
        OC: float
            organic carbon [g/kg]
            
        Returns
        -------
        Bulk density: float    
    """
    
    return 1.449*e**(-0.03*OC)


def hollis_ct(OC, sand, clay):
    """
        Hollis et al., 2012
        n = 333
        RMSE = 0.13, R2 = 0.62

        Parameters
        ----------
        OC: float
            organic carbon [g/kg]

        sand: float 
            Soil volumetric sand content [m**3/m**3]

        clay: float
            Soil volumetric clay content [m**3/m**3]

        Returns
        -------
        Bulk density: float    
    """    
    return 0.80806 + (0.823844*np.exp(-0.27993*OC)) + (0.0014065*sand) - (0.0010299*clay)


def hollis_mh(OC, sand, clay):
    """
        Hollis et al., 2012
        n = 333
        RMSE = 0.13, R2 = 0.62

        Parameters
        ----------
        OC: float
            organic carbon [g/kg]

        sand: float 
            Soil volumetric sand content [m**3/m**3]

        clay: float
            Soil volumetric clay content [m**3/m**3]

        Returns
        -------
        Bulk density: float    
    """    

    return 0.69794 + (0.750636*np.exp(-0.230355*OC)) + (0.0008687*sand) - (0.0005164*clay)


##################### Organic carbon pedotransfer functions ####################################


def orgm2oc(orgm):
    """
        Nelson and Somers, 1996

        Parameters
        ----------
        orgm: float
            Soil volumetric organic matter [%m**3/m**3]
            
        Returns
        -------
        Organic carbon: float  
    """
    return 0.58*orgm


##################### Cation exchange capacity pedotransfer functions ####################################


def krogh(sand, silt, clay, pH_water, OC):
    """
        Krogh et al., 2000

        Parameters
        ----------

        sand: float 
            Soil volumetric sand content [m**3/m**3]

        silt: float
            Soil volumetric silt content [m**3/m**3]

        clay: float
            Soil volumetric clay content [m**3/m**3]

        pH_water: float
            pH of soil water [-]
        
        OC: float
            organic carbon [g/kg]
        
        Returns
        -------
        Cation Exchange Capacity [meq/100g]       
    """
    if (pH_water == np.nan): pH_water = 7.48
    if (OC == np.nan): OC = 0.58
      
    return - 2408.584 + 24.647 * clay + 24.442 * silt + 24.212 * sand - 1.392 * pH_water + 0.243 * OC


def world_wosis_cecph7(silt, clay, OC, depth):
    """
        Mendoza Veirana et al., 2023

        Parameters
        ----------

        clay: float
            Soil volumetric clay content [m**3/m**3]

        silt: float
            Soil volumetric silt content [m**3/m**3]

        OC: float
            organic carbon [g/kg]

        depth: float
            Soil depth [m]
        
        Returns
        -------
        Cation Exchange Capacity [meq/100g]       
    """    

    return -0.07951988335421589 + 4.82413156e-01*clay + 3.59223034e-01*OC + 5.26587906e-02*silt - 8.21084888e-02*depth**2 - 1.01926964e-03*clay**2 -4.24460777e-04*OC**2 -1.74987513e-03*clay*OC -2.09343860e-02*depth*clay + 7.61372784e-02*depth*OC + 2.14532828e-04*depth*silt


################################# # # # # #    PARTICLE DENSITY    # # # # ######################################


def schjonnpd(clay, org, densorg = 1.4, denspart = 2.65, densclay = 2.86):
    """
        Schjonnen et al. 2017

        Parameters
        ----------

        clay: float
            Soil volumetric clay content [m**3/m**3]

        orgm: float
            Soil volumetric organic matter [m**3/m**3]

        Returns
        -------
        Soil particle density [kg/m**3]
    """       
    a = 1.127
    b = 0.373
    c = 2.648
    d = 0.209

    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    return pd

    
##############################################################################################################
############################################# U N C H E C K E D ##############################################
##############################################################################################################

def complete(silt, clay):
    return 100 - silt-clay


def jacob89_pd(clay, org):
    """Schjonnen et al. 2017 using data of Jacobsen 1989. Clay and ORG are expressed in mass fraction"""
    return 2.652 + 0.216 * clay - 2.237 * org

    
def ruhlmann(org, densorg = 1.4, denspart = 2.65, a = 1.127, b = 0.373):
    """
    Rühlmann et al. (2006)
    """
    org = org/100
    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(2.684))**-1
    return pd


def mcbridepd(org):
    """
    McBride et al 2012
    """
    pd = 2.714 - (0.0198*org)
    return pd


###############################################     CLAY CONTENT STIMATIONS     ##########################################


def jacandschj93B_clay(vmc, bd, org, aperm):
    """
    #Jacobsen and Schjonning (1993) Equation (2) for clay content. Here the permittivity is the aparent one.
    """
    return (vmc + 0.0341 - (0.0345)*aperm + (0.00114)*(aperm**2) - (0.0000171)*(aperm**3) + (0.037)*(bd) - (0.00477)*(org))/(0.000736)


###############################################     ORGANIC MATTER STIMATIONS     ##########################################


def devos_inv(bd):
    """
    Best model on De Vos et al. 2005 for soils in Belgium
    """
    return ((-bd+1.775)/0.175)**2


def Banin_Amiel_org(clay, cec):
    """
    CEC [meq/100g] vs clay and organic matter [%], as empirically derived in Banin & amiel 1969 """
    return (-cec + 0.703*clay + 5.1)/2.98


def fernandez88wet(nusell):
    """"Fernandez et al 1988, "Color, Organic Matter, and Pesticide Adsorption Relationships in a Soil Landscape".
    PTF relatinf Nusell values with organic matter [g/kg] for dry and wet samples taken from Indiana, USA. All samples selected had from 58 to 75% silt and from 12 to 35% clay.
    n = 12, r2 = 0.94"""
    return (-nusell+4.38)/0.0523


def fernandez88dry(nusell):
    """"Fernandez et al 1988, "Color, Organic Matter, and Pesticide Adsorption Relationships in a Soil Landscape".
    PTF relatinf Nusell values with organic matter [g/kg] for dry and wet samples taken from Indiana, USA. All samples selected had from 58 to 75% silt and from 12 to 35% clay.
    n = 12, r2 = 0.92"""
    return (-nusell+6.33)/0.0511


######################################## Specific Surface Area (SSA) Pedotransfer ##############################


def Banin_Amiel1(clay):
    """
    Specific Surface Area [m2/g] vs clay content, as empirically derived in Banin & amiel 1969 """
    return 5.76*clay-15.064

""" It is mentioned in Revil etal. 1998 Fgirue 1 an empirical linear relationship among SSA and CEC for rocks"""


def Banin_Amiel2(cec):
    """
    Specific Surface Area [m2/g] vs cec, as empirically derived in Banin & amiel 1969 """
    return (cec-3.23)/0.119


def wang(bound_wat):
    """
    Specific Surface Area [m2/g] vs bound water [% by weight] (called tightly bound water content in this paper) for soils with organic matter. Wang et al. 2011. R2 = 0.827 """
    return 25.2*bound_wat   


    
    ######################################## Bound water Pedotransfer ##############################


def dirksen_dasberg(bd, SSA, l, delta = 3e-10):
    """ 
        Dirksen & Dasberg, 1993
        
        bd: float
            bulk density [g/m3]
            
        SSA: float
            Specific surface area [m2/g]
            
        l: float
            number of molecular layers of bound water
    
    """
    return l*bd*delta*SSA
    
        
######################################## Cation Exchamge Capacity (CEC) Pedotransfer ##############################


def IOPped(clay, org):
    """
    Empirical pedotransfer function fitted to IOP DataSet fro clay and humus vs CEC"""
    return 0.708*clay + 0.762*org - 2.775


def Banin_Amiel3(clay, org):
    """
    CEC [meq/100g] vs clay and organic matter, as empirically derived in Banin & amiel 1969 """
    return 0.703*clay - 2.98*org + 5.1


def shah_singh1(clay):
    """ Shah & Singh empirical model. CEC [meq/100g]"""
    cec = 4.18 + 0.62*clay                            
    return cec
    
    
def bell_keulen95A(clay, org, ph):
    """Bell & Keulen 1995, Soil Pedotransfer Functions for Four Mexican Soils. n = 148 using clay loam (alfisol), sandy loam (entysol), clay (vertisol) and silt loam (alfisol) in Mexico. r2 = 0.94. CEC [cmolc/kg] """
    cec = -10 + 0.163*org*ph - 0.0209*org*clay + 0.131*clay*ph
    return cec
    
    
def bell_keulen95B(clay, org, ph):
    """Bell & Keulen 1995, Soil Pedotransfer Functions for Four Mexican Soils. n = 148 using clay loam (alfisol), sandy loam (entysol), clay (vertisol) and silt loam (alfisol) in Mexico. r2 = 0.96. CEC [cmolc/kg] """
    cec = 42.8 - 5.36*ph +0.297*org - 2.04*clay + 0.363*clay*ph 
    return cec


def mcbratney02(clay, org):
    """McBratney et al. 2002, Eq 10. Empirically developed using a soil database of n=1930 (r2=0.739)aparently exposed in McGarry et al 1989 (McGarry, D., Ward, W.T., McBratney, A.B., 1989. Soil Studies in the Lower Namoi Valley: Methods and Data. The Edgeroi Data Set. CSIRO Division of Soils, Glen Osmond, South Australia.) I did not have accse to. cec[mmol+/kg]. org is actually oragnic carbon."""
    cec = -29.250 + 8.139*clay + 0.253*clay*org
    return cec
    
    
######################################     Cementation exponent (m) Pedotransfer     ########################


def mendelson82(m):
    """Mendelson & Cohen 1982 discrete formula (Eq 28) relating L and m for saturated non-clay rocks 
    with randomly oriented elipsoidal grains for which La + Lb + Lc = 1"""
    l2p = []
    
    for i in range(len(m)):
        l2 = [-3*m[i], 3, -5+3*m[i]]
        roots = np.roots(l2)
        roots = roots[roots.imag == 0 ]
        l2 = roots[roots > 0]
        
        if l2.size==0:
            l2p.append(np.nan)
            
        else:    
            l2p.append(l2[0].real)
            
    return  l2p


def Grunzel(cec):
    """Grunzel 1994 as described in Shah & Singh 2005."""
    return 1.67 + 0.1953*cec**0.5


def schwartza(clay_cont):
    """  F. Schwartza,b,*, Mazdeline E. Schreibera, Tingting Yan. 2008  """
    return 0.485*(clay_cont)**0.0818


################################################################## POROSITY CALCULATION ####################################################################


def park_21(clay, org):
    """
        Park et al, 2021.

        Parameters
        ----------
        clay: float
            clay content [cm3/cm3]

        org: float
            organic matter content [kg/kg]    

        Returns
        -------
        Porosity: float
    """

    return 0.194 + 0.26*clay + 0.65*org


###################################################################    WILTING POINT ##############################################################


def park_19(clay, org):
    """
        Park et al, 2021.

        Parameters
        ----------
        clay: float
            clay content [cm3/cm3]

        org: float
            organic matter content [kg/kg]    

        Returns
        -------
        Wilting point: float
    """

    return 0.02982 + 0.089*clay + 0.65*org
    
    
################################# # # # # QP TO APARENT CONDUCTIVITY # # # # ######################################


def mcneill_eca(qp, frequency, coil_spacing):
    """
    Calculate the apparent electrical conductivity (ECa) based on the raw sensor output (Hs/Hp in ppt).
    >> insert the relevant equation between the square brackets below, 
    >> assigning the output to the variable eca (eca = []) 
    """
    # Constants 
    PI = np.pi
    MU_0 = 4 * PI * 10 ** (-7)
    
    # Angular frequency
    angular_frequency = 2 * PI * frequency

    # Equation to transform the appropriate signal response to ECa, using the low-induction number approximation
    # eca = [] #
    eca = (4 * qp) / (angular_frequency * MU_0 * coil_spacing ** 2)
    return eca


#################################### # # # OUR OWN PTFs PROPOSES  # # # ########################################


def prop1(L):
    a = -1/0.46
    b = 0.71/0.46
    m = ((5/3) - L)/(1 - L**2)
    cc = a/m + b
    return (16.2 + 1.314*cc*100)/1000


def prop2(L):
    a = -1/0.46
    return (((1 - L**2)/(5/3 - L)) - 0.71)*100*a


def prop3(bd, pd, clay_cont, wp, ps = 5, pdn=2.65):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), 
    8, 9 and Archi's law
    And finally we return the value of L
    """
    por = 1- bd/pd
    alpha = -0.46 * clay_cont/100 + 0.71
    m = np.log(((1 - por)*(ps/wp)**alpha) + por) / (np.log(por)*alpha)
    return ((-1/m) + 1)


#def prop5(vmc, clay_cont, organic_cont, _bulk_drydens, pore_ec,wp):
    
#    p = [17.1e-6*(wp**3)*(pore_ec**3), -11.4e-4*(wp**2)*(pore_ec**2), 3.45e-2*(wp)*(pore_ec), 
#         -3.41e-2 - vmc -3.7e-2 * _bulk_drydens + 7.36e-4 * clay_cont + 47.7e-4 * organic_cont]
#    roots = np.roots(p)
#    roots = roots[roots.imag == 0 ]
#    cond = roots[roots > 0]
#    return cond[0].real


def prop6(vmc, por, clay_cont, air_perm, ps, wp):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 4 (for unsaturated consitions), 8, 9 and Archie law for 
    saturated conditions
    And finally we return the value of L
    """
    alpha = -0.46 * clay_cont/100 + 0.71
    
    m = np.log(((1 - por)*(ps/wp)**alpha) + vmc + (por-vmc)*(air_perm/wp)**alpha) / (np.log(por)*alpha)
    
    return ((-1/m) + 1)


def prop7(vmc, por, clay_cont, air_perm, ps, n, wp):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 4 (for unsaturated consitions), 8, 9
    and Archie law for unsaturated conditions.
    Finally we return the value of L
    """
    S     = vmc/(por)
    alpha = -0.46 * clay_cont/100 + 0.71
    
    m     = np.log(((1 - por)*(ps/wp)**alpha) + vmc + (por-vmc)*((air_perm/wp)**alpha) - alpha*n*np.log(S))/ (np.log(por)*alpha)
    
    L     = (-1/m) + 1
    return L
 

def prop8(por, clay_cont, wp, ps = 5):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), 
    8, 9 and Archi's law
    And finally we return the value of L using Han et l. 2020 eq 5
    """
    alpha = -0.46 * clay_cont/100 + 0.71
    m = np.log(((1 - por)*(ps/wp)**alpha) + por) / (np.log(por)*alpha)
    
    p = [-3*m, 3, 3*m-5]
    roots = np.roots(p)
    roots = roots[roots.imag == 0 ]
    cond  = roots[roots > 0]
    
    return cond


def prop9(clay_cont):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), 
    8, 9 and Archi's law
    And finally we return the value of m
    """
    alpha = -0.46 * clay_cont/100 + 0.71
    m     = 1/alpha
    return m


def prop10(bd, pdn, clay_cont, wp, offset):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), 
    8, 9 and Archi's law
    And finally we return the value of m
    """
    rate = 0
    alpha  = -0.46 * clay_cont/100 + 0.71
    por    = 1 - bd/pdn        
    
    m      = (np.log((((1 - por)*(rate)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    return m


def prop11(bd, pdn, clay_cont, perm_sol, wp):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), 
    8, 9 and Archi's law
    And finally we return the value of m
    """
    rate  = perm_sol/wp
    alpha = -0.46 * clay_cont/100 + 0.71
    por   = 1 - bd/pdn    
    
    m     = np.log(((1 - por)*(rate)**alpha) + por) / (np.log(por)*alpha)
    return m


def prop12(bd, pdn, clay_cont, perm_sol, wp, offset):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), Hilhorst 2000 equation and Archi's law
    And finally we return the value of m
    """
    rate   = perm_sol/wp
    alpha  = -0.46 * clay_cont/100 + 0.71
    por    = 1 - bd/pdn                    # deducted from definition of por without assumptions in a three phase medium
    
    m      = (np.log((((1 - por)*(rate)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)
    return m


def prop13(vmc, bd, pdn, clay_cont, perm_sol, wp, offset, n):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (saturated condition for permittivity mixing model), Hilhorst 2000 equation and Archi's law
    And finally we return the value of m
    """
    alpha  = -0.46 * clay_cont/100 + 0.71
    por    = 1 - bd/pdn  
    S      = vmc/por
    rate   = perm_sol/wp
    
    m        = np.log(((((1-por)*np.exp(alpha)) + vmc + ((por-vmc)*(1/wp)**alpha))**(1/alpha)- (offset/wp))*S**(-n))/ np.log(por)
    return m


def prop14(vmc, bd, pdn, clay_cont, perm_sol, wp, m, offset):
    """
    First equation (empirical) relating alpha with clay content is taken from Wunderlich et al. 2013
    The second one is deducted from Brovelli & Casianni 2008 convining eqs. 5 (non-saturated condition for permittivity mixing model), Hilhorst 2000 equation and combined Archi's law. And finally we return the value of n"""
    
    alpha  = -0.46 * clay_cont/100 + 0.71
    por    = 1 - bd/pdn  
    S      = vmc/por
    
    numerator = ((1-por)*perm_sol**alpha + vmc*wp**alpha + (por-vmc))**(1/alpha) 
    
    return np.log((numerator - offset)/((por**m)*wp))/np.log(S)

    
def prop15(clay_cont):
    """Shah & Singh 2005 eq 12b and 12d describing 'm' in funtion of clay content"""
    if clay_cont >= 5:
        m = 0.92*clay_cont**0.2
        
    else:
        m = 1.25
        
    return m


def prop16(clay_cont):
    """
    Combination of Grunzet thesis 1994 pedotransfer function (I could not access to) and Shah & Singh 2005 """
    cec = 4.18 + 0.62*clay_cont                # Fig 5, Shah & Singh 2005 [meq/100g]
    m   = 1.6 + 0.1953*cec**0.5                # Gunzel 1994 as described in Shah & Singh 2005. Original expression: 1.6+0.1953*cec**0.5  
    return m

def prop17(bd, part_dens, ps, wp, CEC, offset):
    """
     Eq. 16"""

    por = 1 - bd/part_dens
    alpha =-0.00076811*CEC**2 + 0.05330303*CEC + 0.40328599
    m = (np.log((((1 - por)*(ps/wp)**alpha) + por)**(1/alpha) - (offset/wp))) / np.log(por)               
    
    return m