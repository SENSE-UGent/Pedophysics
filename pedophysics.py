#################
import pedotransfer_functions as ptf
import pedophysical_models as pm
import pandas as pd
import numpy as np

#################
class Soil(object):
    """ A virtual soil that can be composed of several properties """
    def __init__(
            self, n_states = np.nan, temperature = np.nan, water = np.nan, bound_water = np.nan, water_init = np.nan, salinity = np.nan, 
            sand = np.nan, silt = np.nan, clay = np.nan, bulk_density = np.nan, particle_density = np.nan, CEC = np.nan, LOI = np.nan, OC = np.nan, 
            orgm = np.nan, 
            bulk_perm = np.nan, bulk_iperm = np.nan, bulk_aperm = np.nan, bulk_perm_inf = np.nan, bulk_perm_init = np.nan, 
            water_perm = np.nan, bound_water_perm = np.nan, water_iperm = np.nan, bound_water_iperm = np.nan, 
            solid_perm = np.nan, air_perm = np.nan, offset = np.nan, 
            bulk_cond = np.nan, bulk_icond = np.nan, water_cond = np.nan, matrix_cond = np.nan, air_cond = np.nan,
            L = np.nan, Lw = np.nan, m = np.nan, n = np.nan, alpha = np.nan, 
            frequency_perm = np.nan, frequency_cond = np.nan,
            land = np.nan, region = np.nan, depth = np.nan, texture = np.nan, instrument = np.nan, pH_water = np.nan, df = np.nan):

        """
            Pedophysics
            Soil object
        
            Parameters
            ----------
            temperature: float
                Soil temperature in celcius [C]

            water: float
                Soil volumetric water content [m**3/m**3]

            bound_water: float
                Soil volumetric bound water content [m**3/m**3] 
            
            water_init: float
                Soil minimum volumetric water content between all the states [m**3/m**3] 
            
            salinity: float
                Soil (NaCl) salinity of the bulk pore fluid [mol/L]
        
            sand: float
                Soil volumetric sand content [m**3/m**3]
        
            silt:  float
                Soil volumetric silt content [m**3/m**3]
        
            clay: float
                Soil volumetric clay content [m**3/m**3]
        
            bulk_density: float
                Soil bulk density [kg/m**3]
        
            particle_density: float
                Soil particle density [kg/m**3]
        
            CEC: float
                Soil cation exchange capacity [meq/100g]

            LOI: float
                Soil loss on ignition [g]?
        
            OC: float
                Soil organic carbon [%g/g] ?
        
            orgm: float
                Soil volumetric organic matter [%m**3/m**3]

            pH_water: float
                pH of the soil water

            bulk_perm: float
                Soil bulk real relative dielectric permittivity [-]
        
            bulk_iperm: float
                Soil bulk imaginary relative dielectric permittivity [-]

            bulk_aperm: float
                Soil bulk apparent relative dielectric permittivity [-]
        
            bulk_perm_inf: float
                Soil bulk real relative permittivity at infinite frequency [-]
                
            bulk_perm_init: float
                Soil minimun bulk real relative permittivity [-]

            water_perm: float
                Soil water phase real dielectric permittivity [-]
        
            bound_water_perm: float
                Soil bound water phase real relative dielectric permittivity [-]
        
            water_iperm: float
                Soil water phase imaginary relative dielectric permittivity [-]

            bound_water_iperm: float
                Soil bound water phase imaginary relative dielectric permittivity [-]
        
            solid_perm: float
                Soil solid real dielectric permittivity phase [-]
        
            air_perm: float
                Soil air real dielectric permittivity phase [-]
        
            offset: float
                Soil offset as defined in Hilhorst (2000) formula [-]
        
            bulk_cond: float
                Soil bulk real electrical conductivity [S/m]
        
            bulk_icond: float
                Soil bulk imaginary electrical conductivity [S/m]
        
            water_cond: float
                Soil water real electrical conductivity [S/m]
        
            matrix_cond: float
                Soil maxtrix real electrical conductivity [S/m]
        
            air_cond: float
                Soil air real electrical conductivity [S/m]
        
            L: float
                Soil particle depolarization factor as defined in effective medium theory [-]
        
            Lw: float
                Soil water depolarization factor as defined in effective medium theory [-]
        
            m: float
                Soil cementation factor as defined in Archie law [-]

            n: float
                Soil saturation factor as defined in Archie second law [-]
        
            alpha: float
                Soil alpha exponent as defined in Roth et al. (1990)
        
            frequency_perm: float
                Frequency of dielectric permittivity measurement [Hz]
        
            frequency_cond: float
                Frequency of electric conductivity measurement [Hz]
        
            land: str
                Land use: 'forest', 'pasture', 'cultivated topsoil', 'desert', 'archaeology'

            region : str
                Region of the world: 'Europe', "North America", "Central America", "South America"

            depth: float
                Soil depth [m] 

            texture: str
                Soil texture according to USDA convention: "Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Silty clay loam", "Sandy clay", "Clay", "Silty clay"  

            instrument: str
                Intrument utilized: 'TDR', 'GPR', 'HydraProbe'
            
                
            Returns
            -------
            Soil: self object    
        """

        assert (type(temperature) == float) or (type(temperature) == np.ndarray) or (type(temperature) == int) or (temperature == np.nan), "'temperature' must be a int, float, numpy.array or np.nan"
        assert (type(water) == float) or (type(water) == np.ndarray) or (type(water) == int) or (water == np.nan), "'water' must be an int, float, numpy.array or np.nan"
        assert (type(bound_water) == float) or (type(bound_water) == np.ndarray) or (bound_water == np.nan), "'bound_water' must be a float, numpy.array or np.nan"
        assert (type(water_init) == float) or (type(water_init) == int) or (water_init == np.nan), "'water_init' must be an int, float or np.nan"
        assert (type(salinity) == float) or (type(salinity) == np.ndarray) or (type(salinity) == int) or (salinity == np.nan), "'salinity' must be an int, float, numpy.array or np.nan"
        assert (type(sand) == float) or (type(sand) == int) or (type(sand) == np.ndarray) or(sand == np.nan), "'sand' must be an int, float or np.nan"
        assert (type(silt) == float) or (type(silt) == int) or (type(silt) == np.ndarray) or (silt == np.nan), "'silt' must be an int, float or np.nan"
        assert (type(clay) == float) or (type(clay) == int) or (type(clay) == np.ndarray) or (clay == np.nan), "'clay' must be an int, float or np.nan"
        assert (type(bulk_density) == float) or (type(bulk_density) == np.ndarray) or (type(bulk_density) == int) or (bulk_density == np.nan), "'bulk_density' must be an int, float, numpy.array or np.nan"
        assert (type(particle_density) == float) or (type(particle_density) == np.ndarray) or (type(particle_density) == int) or (type(particle_density) == list) or (particle_density == np.nan), "'particle_density' must be an int, float, numpy.array, list or np.nan"
        assert (type(CEC) == float) or (type(CEC) == np.ndarray) or (type(CEC) == int) or (type(CEC) == list) or (CEC == np.nan), "'CEC' must be an int, float, numpy.array, list or np.nan"
        assert (type(LOI) == float) or (type(LOI) == np.ndarray) or (type(LOI) == int)  or (type(LOI) == list) or (LOI == np.nan), "'LOI' must be an int, float, numpy.array, list or np.nan"
        assert (type(OC) == float) or (type(OC) == np.ndarray) or (type(OC) == int) or (type(OC) == list) or (OC == np.nan), "'OC' must be an int, float, numpy.array, list or np.nan"
        assert (type(orgm) == float) or (type(orgm) == np.ndarray) or (type(orgm) == int) or (type(orgm) == list) or (orgm == np.nan), "'orgm' must be an int, float, numpy.array, list or np.nan"
        assert (type(pH_water) == float) or (type(pH_water) == np.ndarray) or (type(pH_water) == int) or (type(pH_water) == list) or (pH_water == np.nan), "'pH_water' must be an int, float, numpy.array, list or np.nan"
        assert (type(bulk_perm) == float) or (type(bulk_perm) == np.ndarray) or (type(bulk_perm) == int) or (bulk_perm == np.nan), "'bulk_perm' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_iperm) == float) or (type(bulk_iperm) == np.ndarray) or (type(bulk_iperm) == int) or (bulk_iperm == np.nan), "'bulk_iperm' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_aperm) == float) or (type(bulk_aperm) == np.ndarray) or (type(bulk_aperm) == int) or (bulk_aperm == np.nan), "'bulk_aperm' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_perm_inf) == float) or (type(bulk_perm_inf) == np.ndarray) or (type(bulk_perm_inf) == int) or (bulk_perm_inf == np.nan), "'bulk_perm_inf' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_perm_init) == float) or (type(bulk_perm_init) == int) or (bulk_perm_init == np.nan), "'bulk_perm_init' must be an int, float or np.nan"
        assert (type(water_perm) == float) or (type(water_perm) == np.ndarray) or (type(water_perm) == int) or (water_perm == np.nan), "'water_perm' must be an int, float, numpy.array or np.nan"
        assert (type(bound_water_perm) == float) or (type(bound_water_perm) == np.ndarray) or (type(bound_water_perm) == int) or (bound_water_perm == np.nan), "'bound_water_perm' must be an int, float, numpy.array or np.nan"
        assert (type(water_iperm) == float) or (type(water_iperm) == np.ndarray) or (type(water_iperm) == int) or (water_iperm == np.nan), "'water_iperm' must be an int, float, numpy.array or np.nan"
        assert (type(bound_water_iperm) == float) or (type(bound_water_iperm) == np.ndarray) or (type(bound_water_iperm) == int) or (bound_water_iperm == np.nan), "'bound_water_iperm' must be an int, float, numpy.array or np.nan"
        assert (type(solid_perm) == float) or (type(solid_perm) == np.ndarray) or (type(solid_perm) == int) or (solid_perm == np.nan), "'solid_perm' must be an int, float, numpy.array or np.nan"
        assert (type(air_perm) == float) or (type(air_perm) == np.ndarray) or (type(air_perm) == int) or (air_perm == np.nan), "'air_perm' must be an int, float, numpy.array or np.nan"
        assert (type(offset) == float) or (type(offset) == np.ndarray) or (type(offset) == int) or (offset == np.nan), "'offset' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_cond) == float) or (type(bulk_cond) == np.ndarray) or (type(bulk_cond) == int) or (bulk_cond == np.nan), "'bulk_cond' must be an int, float, numpy.array or np.nan"
        assert (type(bulk_icond) == float) or (type(bulk_icond) == np.ndarray) or (type(bulk_icond) == int) or (bulk_icond == np.nan), "'bulk_icond' must be an int, float, numpy.array or np.nan"
        assert (type(water_cond) == float) or (type(water_cond) == np.ndarray) or (type(water_cond) == int) or (water_cond == np.nan), "'water_cond' must be an int, float, numpy.array or np.nan"
        assert (type(matrix_cond) == float) or (type(matrix_cond) == np.ndarray) or (type(matrix_cond) == int) or (matrix_cond == np.nan), "'matrix_cond' must be an int, float, numpy.array or np.nan"
        assert (type(air_cond) == float) or (type(air_cond) == np.ndarray) or (type(air_cond) == int) or (air_cond == np.nan), "'air_cond' must be an int, float, numpy.array or np.nan"
        assert (type(L) == float) or (type(L) == np.ndarray) or (type(L) == int) or (L == np.nan), "'L' must be an int, float, numpy.array or np.nan"
        assert (type(Lw) == float) or (type(Lw) == np.ndarray) or (type(Lw) == int) or (Lw == np.nan), "'Lw' must be an int, float, numpy.array or np.nan"
        assert (type(m) == float) or (type(m) == np.ndarray) or (type(m) == int) or (m == np.nan), "'m' must be an int, float, numpy.array or np.nan"
        assert (type(n) == float) or (type(n) == np.ndarray) or (type(n) == int) or (n == np.nan), "'n' must be an int, float, numpy.array or np.nan"
        assert (type(alpha) == float) or (type(alpha) == np.ndarray) or (type(alpha) == int) or (alpha == np.nan), "'alpha' must be an int, float, numpy.array or np.nan"
        assert (type(frequency_perm) == float) or (type(frequency_perm) == np.ndarray) or (type(frequency_perm) == int) or (frequency_perm == np.nan), "'frequency_perm' must be an int, float, numpy.array or np.nan"
        assert (type(frequency_cond) == float) or (type(frequency_cond) == np.ndarray) or (type(frequency_cond) == int) or (frequency_cond == np.nan), "'frequency_cond' must be an int, float, numpy.array or np.nan"
        #assert (type(land) == str) or (land == np.nan), "'land' must be a str or np.nan"
        #assert (type(region) == str or (region == np.nan)), "'region' must be a str or np.nan"
        assert (type(depth) == float) or (type(depth) == int) or (depth == np.nan), "'depth' must be an int, float or np.nan. 'A soil cannot be at the same time in two different places'"
        assert texture in ["Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Sandy clay", "Clay", "Silty clay", np.nan], "texture must be 'Sand' or 'Loamy sand' or 'Sandy loam' or 'Loam' or 'Silt loam' or 'Silt' or 'Sandy clay loam' or 'Clay loam' or 'Sandy clay' or 'Clay' or 'Silty clay' or np.nan"
        assert instrument in ["TDR", "GPR", 'HydraProbe', np.nan], "'instrument' must be 'TDR', 'GPR', 'HydraProbe' or np.nan"


### Here we define the state variables
        self.temperature = temperature                      # Soil temperature in celcius [C]
        if type(self.temperature) != np.ndarray:
            self.temperature = np.array([temperature])  

        self.water = water                                  # Soil volumetric water content [m**3/m**3]
        if type(self.water) != np.ndarray:
            self.water = np.array([water])    

        self.bound_water = bound_water                      # Soil volumetric bound water content [m**3/m**3] 
        if type(self.bound_water) != np.ndarray:
            self.bound_water = np.array([bound_water])  

        self.salinity = salinity                            # Soil (NaCl) salinity of the bulk pore fluid [mol/L]
        if type(self.salinity) != np.ndarray:
            self.salinity = np.array([salinity])  

        self.sand = sand                                    # Soil volumetric sand content [m**3/m**3]
        if type(self.sand) != np.ndarray:
            self.sand = np.array([sand]) 

        self.silt = silt                                    # Soil volumetric silt content [m**3/m**3]
        if type(self.silt) != np.ndarray:
            self.silt = np.array([silt])

        self.clay = clay                                    # Soil volumetric clay content [m**3/m**3]
        if type(self.clay) != np.ndarray:
            self.clay = np.array([clay])

        self.bulk_density = bulk_density                    # Soil bulk density [g/cm3]
        if type(self.bulk_density) != np.ndarray:
            self.bulk_density = np.array([bulk_density])  

        self.particle_density = particle_density            # Soil particle density [g/cm3]
        if type(self.particle_density) == list:
            self.particle_density = np.array(particle_density)  
        if type(self.particle_density) != np.ndarray:
            self.particle_density = np.array([particle_density])  

        self.CEC = CEC                                      # Soil cation exchange capacity [meq/100g]
        if type(self.CEC) == list:
            self.CEC = np.array(CEC)
        if type(self.CEC) != np.ndarray:
            self.CEC = np.array([CEC])  

        self.LOI = LOI                                      # Soil loss on ignition [g]?
        if type(self.LOI) == list:
            self.LOI = np.array(LOI) 
        if type(self.LOI) != np.ndarray:
            self.LOI = np.array([LOI]) 

        self.OC = OC                                        # Soil organic carbon [%g/g] ?
        if type(self.OC) == list:
            self.OC = np.array(OC)
        if type(self.OC) != np.ndarray:
            self.OC = np.array([OC])

        self.orgm = orgm                                    # Soil volumetric organic matter [%m**3/m**3]
        if type(self.orgm) == list:
            self.orgm = np.array(orgm)
        if type(self.orgm) != np.ndarray:
            self.orgm = np.array([orgm])

        self.pH_water = pH_water                            # pH of the soil water
        if type(self.pH_water) == list:
            self.pH_water = np.array(pH_water)
        if type(self.pH_water) != np.ndarray:
            self.pH_water = np.array([pH_water])

        self.bulk_perm = bulk_perm                          # Soil bulk real relative dielectric permittivity [-]
        if type(self.bulk_perm) != np.ndarray:
            self.bulk_perm = np.array([bulk_perm])      

        self.bulk_iperm = bulk_iperm                        # Soil bulk imaginary relative dielectric permittivity [-]
        if type(self.bulk_iperm) != np.ndarray:
            self.bulk_iperm = np.array([bulk_iperm])  

        self.bulk_aperm = bulk_aperm                        # Soil bulk apparent relative dielectric permittivity [-]
        if type(self.bulk_aperm) != np.ndarray:
            self.bulk_aperm = np.array([bulk_aperm])  

        self.bulk_perm_inf = bulk_perm_inf                  # Soil bulk real relative dielectric permittivity at infinite frequency [-]
        if type(self.bulk_perm_inf) != np.ndarray:
            self.bulk_perm_inf = np.array([bulk_perm_inf])  

        self.water_perm = water_perm                        # Soil water phase real relative dielectric permittivity [-]
        if type(self.water_perm) != np.ndarray:
            self.water_perm = np.array([water_perm])  

        self.water_iperm = water_iperm                      # Soil water phase imaginary relative dielectric permittivity [-]
        if type(self.water_iperm) != np.ndarray:
            self.water_iperm = np.array([water_iperm])  

        self.bulk_cond = bulk_cond                          # Soil bulk real electrical conductivity [S/m]
        if type(self.bulk_cond) != np.ndarray:
            self.bulk_cond = np.array([bulk_cond])  

        self.bulk_icond = bulk_icond                        # Soil bulk imaginary electrical conductivity [S/m]
        if type(self.bulk_icond) != np.ndarray:
            self.bulk_icond = np.array([bulk_icond])  

        self.water_cond = water_cond                        # Soil water real electrical conductivity [S/m]
        if type(self.water_cond) != np.ndarray:
            self.water_cond = np.array([water_cond])  

        self.frequency_perm = frequency_perm                # Frequency of soil dielectric permittivity measurement [Hz]
        if type(self.frequency_perm) != np.ndarray:
            self.frequency_perm = np.array([frequency_perm])
            
        self.frequency_cond = frequency_cond                # Frequency of soil electric conductivity measurement [Hz]
        if type(self.frequency_cond) != np.ndarray:
            self.frequency_cond = np.array([frequency_cond])


### Here we define non-state variables
        self.water_init = water_init                        # Soil minimum volumetric water content between all the states [m**3/m**3] 
        self.bulk_perm_init = bulk_perm_init                # Soil minimun bulk real relative permittivity
        self.bound_water_perm = bound_water_perm            # Soil bound water phase real relative dielectric permittivity [-]
        self.bound_water_iperm = bound_water_iperm          # Soil bound water phase imaginary relative dielectric permittivity [-]
        self.solid_perm = solid_perm                        # Soil solid real dielectric permittivity phase [-]
        self.air_perm = air_perm                            # Soil air real dielectric permittivity phase [-]
        self.offset = offset                                # Soil offset as defined in Hilhorst (2000) formula [-]
        self.matrix_cond = matrix_cond                      # Soil maxtrix real electrical conductivity [S/m]
        self.air_cond = air_cond                            # Soil air real electrical conductivity [S/m]
        self.L = L                                          # Soil particle depolarization factor as defined in effective medium theory [-]
        self.Lw = Lw                                        # Soil water depolarization factor as defined in effective medium theory [-]
        self.m = m                                          # Soil cementation factor as defined in Archie law [-]
        self.n = n                                          # Soil saturation factor as defined in Archie second law [-]
        self.alpha = alpha                                  # Soil alpha exponent as defined in Roth et al. (1990)       
        self.land = land                                    # Land use: 'forest', 'pasture', 'cultivated topsoil', 'desert', 'archaeology'
        self.region = region                                # Region of the world: 'Europe', "North America", "Central America", "South America"
        self.depth = depth                                  # Soil depth [m] 
        self.texture = texture                              # Soil texture: "Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Sandy clay", "Clay", "Silty clay" 
        self.instrument = instrument                        # Intrument utilized: 'TDR', 'GPR', 'HydraProbe'

        n_states = max(   len(self.temperature), len(self.water), len(self.bound_water), len(self.salinity), 
        len(self.sand), len(self.silt), len(self.clay), len(self.bulk_density), len(self.particle_density),
        len(self.CEC), len(self.LOI), len(self.OC), len(self.orgm), len(self.pH_water),  
        len(self.bulk_perm), len(self.bulk_iperm), len(self.bulk_aperm), len(self.bulk_perm_inf), len(self.water_perm), 
        len(self.water_iperm), len(self.bulk_cond), len(self.bulk_icond), len(self.water_cond), 
        len(self.frequency_perm), len(self.frequency_cond) )

        self.n_states = n_states                            # Number of states of the soil


### Fill the state variables with nans when are shorter than n_states
        if len(self.temperature) != n_states :
            self.temperature = np.append(self.temperature, [np.nan]*(n_states - len(self.temperature)))
        
        if len(self.water) != n_states:
            self.water = np.append(self.water, [np.nan]*(n_states - len(self.water)))    

        if len(self.bound_water) != n_states:
            self.bound_water = np.append(self.bound_water, [np.nan]*(n_states - len(self.bound_water)))   

        if len(self.salinity) != n_states:
            self.salinity = np.append(self.salinity, [np.nan]*(n_states - len(self.salinity)))   

        if len(self.sand) != n_states:
            self.sand = np.append(self.sand, [np.nan]*(n_states - len(self.sand))) 

        if len(self.silt) != n_states:
            self.silt = np.append(self.silt, [np.nan]*(n_states - len(self.silt))) 

        if len(self.clay) != n_states:
            self.clay = np.append(self.clay, [np.nan]*(n_states - len(self.clay))) 

        if len(self.bulk_density) != n_states:
            self.bulk_density = np.append(self.bulk_density, [np.nan]*(n_states - len(self.bulk_density)))  

        if len(self.particle_density) != n_states:
            self.particle_density = np.append(self.particle_density, [np.nan]*(n_states - len(self.particle_density))) 

        if len(self.CEC) != n_states:
            self.CEC = np.append(self.CEC, [np.nan]*(n_states - len(self.CEC)))

        if len(self.LOI) != n_states:
            self.LOI = np.append(self.LOI, [np.nan]*(n_states - len(self.LOI)))  

        if len(self.OC) != n_states:
            self.OC = np.append(self.OC, [np.nan]*(n_states - len(self.OC)))

        if len(self.orgm) != n_states:
            self.orgm = np.append(self.orgm, [np.nan]*(n_states - len(self.orgm)))    

        if len(self.bulk_perm ) != n_states:
            self.bulk_perm  = np.append(self.bulk_perm , [np.nan]*(n_states - len(self.bulk_perm )))  
                
        if len(self.bulk_iperm ) != n_states:
            self.bulk_iperm  = np.append(self.bulk_iperm , [np.nan]*(n_states - len(self.bulk_iperm )))  

        if len(self.bulk_aperm ) != n_states:
            self.bulk_aperm  = np.append(self.bulk_aperm , [np.nan]*(n_states - len(self.bulk_aperm )))  

        if len(self.bulk_perm_inf ) != n_states:
            self.bulk_perm_inf  = np.append(self.bulk_perm_inf , [np.nan]*(n_states - len(self.bulk_perm_inf )))  

        if len(self.water_perm ) != n_states:
            self.water_perm  = np.append(self.water_perm , [np.nan]*(n_states - len(self.water_perm )))  

        if len(self.water_iperm ) != n_states:
            self.water_iperm  = np.append(self.water_iperm , [np.nan]*(n_states - len(self.water_iperm )))  
 
        if len(self.bulk_cond  ) != n_states:
            self.bulk_cond   = np.append(self.bulk_cond  , [np.nan]*(n_states - len(self.bulk_cond  )))  

        if len(self.bulk_icond  ) != n_states:
            self.bulk_icond   = np.append(self.bulk_icond  , [np.nan]*(n_states - len(self.bulk_icond  )))  

        if len(self.water_cond  ) != n_states:
            self.water_cond   = np.append(self.water_cond  , [np.nan]*(n_states - len(self.water_cond ))) 
 
        if len(self.frequency_perm  ) != n_states:
            self.frequency_perm   = np.append(self.frequency_perm, [np.nan]*(n_states - len(self.frequency_perm ))) 

        if len(self.frequency_cond ) != n_states:
            self.frequency_cond   = np.append(self.frequency_cond, [np.nan]*(n_states - len(self.frequency_cond ))) 


### If just the state zero is defined in a state variable, then fully fill it.
        if (self.temperature[0] != np.nan) and (self.temperature[1:(n_states)] == np.nan).all():
            self.temperature[1:n_states] = self.temperature[0] 

        if (self.water[0] != np.nan) and (self.water[1:(n_states)] == np.nan).all():
            self.water[1:n_states] = self.water[0] 

        if (self.bound_water[0] != np.nan) and (self.bound_water[1:(n_states)] == np.nan).all():
            self.bound_water[1:n_states] = self.bound_water[0] 

        if (self.salinity[0] != np.nan) and (self.salinity[1:(n_states)] == np.nan).all():
            self.salinity[1:n_states] = self.salinity[0] 

        if (self.sand[0] != np.nan) and (self.sand[1:(n_states)] == np.nan).all():
            self.sand[1:n_states] = self.sand[0] 

        if (self.silt[0] != np.nan) and (self.silt[1:(n_states)] == np.nan).all():
            self.silt[1:n_states] = self.silt[0] 

        if (self.clay[0] != np.nan) and (self.clay[1:(n_states)] == np.nan).all():
            self.clay[1:n_states] = self.clay[0] 

        if (self.bulk_density[0] != np.nan) and (np.isnan(self.bulk_density[1:(n_states)])).all():
            self.bulk_density[1:n_states] = self.bulk_density[0] 

        if (self.particle_density[0] != np.nan) and (np.isnan(self.particle_density[1:(n_states)])).all():
            self.particle_density[1:n_states] = self.particle_density[0] 

        if (self.CEC[0] != np.nan) and (self.CEC[1:(n_states)] == np.nan).all():
            self.CEC[1:n_states] = self.CEC[0] 

        if (self.LOI[0] != np.nan) and (self.LOI[1:(n_states)] == np.nan).all():
            self.LOI[1:n_states] = self.LOI[0] 

        if (self.OC[0] != np.nan) and (self.OC[1:(n_states)] == np.nan).all():
            self.OC[1:n_states] = self.OC[0] 

        if (self.orgm[0] != np.nan) and (self.orgm[1:(n_states)] == np.nan).all():
            self.orgm[1:n_states] = self.orgm[0] 

        if (self.bulk_perm[0] != np.nan) and (np.isnan(self.bulk_perm[1:(n_states)])).all():
            self.bulk_perm[1:n_states] = self.bulk_perm[0] 

        if (self.bulk_iperm[0] != np.nan) and (np.isnan(self.bulk_iperm[1:(n_states)])).all():
            self.bulk_iperm[1:n_states] = self.bulk_iperm[0] 

        if (self.bulk_aperm[0] != np.nan) and (np.isnan(self.bulk_aperm[1:(n_states)])).all():
            self.bulk_aperm[1:n_states] = self.bulk_aperm[0] 

        if (self.bulk_perm_inf[0] != np.nan) and (np.isnan(self.bulk_perm_inf[1:(n_states)])).all():
            self.bulk_perm_inf[1:n_states] = self.bulk_perm_inf[0] 

        if (self.water_perm[0] != np.nan) and (np.isnan(self.water_perm[1:(n_states)])).all():
            self.water_perm[1:n_states] = self.water_perm[0] 

        if (self.water_iperm[0] != np.nan) and (np.isnan(self.water_iperm[1:(n_states)])).all():
            self.water_iperm[1:n_states] = self.water_iperm[0] 

        if (self.bulk_cond[0] != np.nan) and (np.isnan(self.bulk_cond[1:(n_states)])).all():
            self.bulk_cond[1:n_states] = self.bulk_cond[0] 

        if (self.bulk_icond[0] != np.nan) and (np.isnan(self.bulk_icond[1:(n_states)])).all():
            self.bulk_icond[1:n_states] = self.bulk_icond[0] 

        if (self.water_cond[0] != np.nan) and (np.isnan(self.water_cond[1:(n_states)])).all():
            self.water_cond[1:n_states] = self.water_cond[0] 

        if (self.frequency_perm[0] != np.nan) and (np.isnan(self.frequency_perm[1:(n_states)])).all():
            self.frequency_perm[1:n_states] = self.frequency_perm[0] 

        if (self.frequency_cond[0] != np.nan) and (np.isnan(self.frequency_cond[1:(n_states)])).all():
            self.frequency_cond[1:n_states] = self.frequency_cond[0] 


### Defining Soil.df  
        df = pd.DataFrame({
            'temperature': self.temperature,
            'water': self.water,
            'bound_water': self.bound_water,
            'salinity': self.salinity,
            'sand': self.sand,
            'silt': self.silt,
            'clay': self.clay,
            'bulk_density': self.bulk_density,
            'particle_density': self.particle_density,
            'CEC': self.CEC,
            'LOI': self.LOI,
            'OC': self.OC,
            'orgm': self.orgm,
            'pH_water' : pH_water,
            'bulk_perm': self.bulk_perm,
            'bulk_iperm': self.bulk_iperm,
            'bulk_aperm': self.bulk_aperm,
            'bulk_perm_inf': self.bulk_perm_inf,
            'water_perm': self.water_perm,
            'water_iperm': self.water_iperm,
            'bulk_cond': self.bulk_cond,
            'bulk_icond': self.bulk_icond,
            'water_cond': self.water_cond,
            'frequency_perm': self.frequency_perm,
            'frequency_cond': self.frequency_cond })

        self.df = df

    def get_n_states(self):
        return self.n_states
    def get_temperature(self):
        return self.temperature
    def get_water(self):
        return self.water
    def get_bound_water(self):
        return self.bound_water
    def get_water_init(self):
        return self.water_init
    def get_salinity(self):
        return self.salinity
    def get_sand(self):
        return self.sand
    def get_silt(self):
        return self.silt
    def get_clay(self):
        return self.clay
    def get_bulk_density(self):
        return self.bulk_density
    def get_particle_density(self):
        return self.particle_density
    def get_CEC(self):
        return self.CEC
    def get_LOI(self):
        return self.LOI
    def get_OC(self):
        return self.OC
    def get_orgm(self):
        return self.orgm
    def get_bulk_perm(self):
        return self.bulk_perm
    def get_bulk_iperm(self):
        return self.bulk_iperm
    def get_bulk_aperm(self):
        return self.bulk_aperm
    def get_bulk_perm_inf(self):
        return self.bulk_perm_inf
    def get_bulk_perm_init(self):
        return self.bulk_perm_init
    def get_water_perm(self):
        return self.water_perm
    def get_bound_water_perm(self):
        return self.bound_water_perm
    def get_water_iperm(self):
        return self.water_iperm
    def get_bound_water_iperm(self):
        return self.bound_water_iperm
    def get_solid_perm(self):
        return self.solid_perm
    def get_air_perm(self):
        return self.air_perm
    def get_offset(self):
        return self.offset
    def get_bulk_cond(self):
        return self.bulk_cond
    def get_bulk_icond(self):
        return self.bulk_icond
    def get_water_cond(self):
        return self.water_cond
    def get_matrix_cond(self):
        return self.matrix_cond
    def get_air_cond(self):
        return self.air_cond
    def get_L(self):
        return self.L
    def get_Lw(self):
        return self.Lw
    def get_m(self):
        return self.m
    def get_n(self):
        return self.n
    def get_alpha(self):
        return self.alpha
    def get_frequency_perm(self):
        return self.frequency_perm
    def get_frequency_cond(self):
        return self.frequency_cond
    def get_land(self):
        return self.land
    def get_region(self):
        return self.region
    def get_depth(self):
        return self.depth   
    def get_texture(self):
        return self.texture 
    def get_instrument(self):
        return self.instrument 
    def get_pH_water(self):
        return self.pH_water 

    def set_n_states(self, new_n_states):
        self.n_states = new_n_states
    def set_temperature(self, new_temperature):
        self.temperature= new_temperature
    def set_water(self, new_water):
        self.water = new_water
    def set_bound_water(self, new_bound_water):
        self.bound_water = new_bound_water
    def set_water_init(self, new_water_init):
        self.water_init = new_water_init 
    def set_salinity(self, new_salinity):
        self.salinity = new_salinity
    def set_sand(self, new_sand):
        self.sand = new_sand
    def set_silt(self, new_silt):
        self.silt = new_silt
    def set_clay(self, new_clay):
        self.clay = new_clay
    def set_bulk_density(self, new_bulk_density):
        self.bulk_density = new_bulk_density
    def set_particle_density(self, new_particle_density):
        self.particle_density = new_particle_density
    def set_CEC(self, new_CEC):
        self.CEC = new_CEC
    def set_LOI(self, new_LOI):
        self.LOI = new_LOI
    def set_OC(self, new_OC):
        self.OC = new_OC
    def set_orgm(self, new_orgm):
        self.orgm = new_orgm
    def set_bulk_perm(self, new_bulk_perm):
        self.bulk_perm = new_bulk_perm
    def set_bulk_iperm(self, new_bulk_iperm):
        self.bulk_iperm = new_bulk_iperm
    def set_bulk_aperm(self, new_bulk_aperm):
        self.bulk_aperm = new_bulk_aperm
    def set_bulk_perm_inf(self, new_bulk_perm_inf):
        self.bulk_perm_inf = new_bulk_perm_inf
    def set_bulk_perm_init(self, new_bulk_perm_init):
        self.bulk_perm_init = new_bulk_perm_init
    def set_water_perm(self, new_water_perm):
        self.water_perm = new_water_perm
    def set_bound_water_perm(self, new_bound_water_perm):
        self.bound_water_perm = new_bound_water_perm
    def set_water_iperm(self, new_water_iperm):
        self.water_iperm = new_water_iperm
    def set_bound_water_iperm(self, new_bound_water_iperm):
        self.bound_water_iperm = new_bound_water_iperm
    def set_solid_perm(self, new_solid_perm):
        self.solid_perm = new_solid_perm
    def set_air_perm(self, new_air_perm):
        self.air_perm = new_air_perm
    def set_offset(self, new_offset):
        self.offset = new_offset
    def set_bulk_cond(self, new_bulk_cond):
        self.bulk_cond = new_bulk_cond
    def set_bulk_icond(self, new_bulk_icond):
        self.bulk_icond = new_bulk_icond
    def set_water_cond(self, new_water_cond):
        self.water_cond = new_water_cond
    def set_matrix_cond(self, new_matrix_cond):
        self.matrix_cond = new_matrix_cond
    def set_air_cond(self, new_air_cond):
        self.air_cond = new_air_cond
    def set_L(self, new_L):
        self.L = new_L
    def set_Lw(self, new_Lw):
        self.Lw = new_Lw
    def set_m(self, new_m):
        self.m = new_m
    def set_n(self, new_n):
        self.n = new_n
    def set_alpha(self, new_alpha):
        self.alpha = new_alpha
    def set_frequency_perm(self, new_frequency_perm):
        self.frequency_perm = new_frequency_perm
    def set_frequency_cond(self, new_frequency_cond):
        self.frequency_cond = new_frequency_cond
    def set_land(self, new_land):
        self.land = new_land
    def set_region(self, new_region):
        self.region = new_region
    def set_depth(self, new_depth):
        self.depth = new_depth
    def set_texture(self, new_texture):
        self.texture = new_texture
    def set_instrument(self, new_instrument):
        self.instrument = new_instrument
    def set_pH_water(self, new_pH_water):
        self.pH_water = new_pH_water
        
        
    def __str__(self):                                                   
        """ Returns a string representation of self """
        return str(self.df)



####  TODO function 'save' to export a .json

# "<<<" + 'n_states =' + str(self.n_states) + '[-]' + ", Temperature = " + str(self.temperature) + " Celsius" + ", Water content = " + str(self.water) + " [m**3/m**3]" + ", Bound water = " + str(self.bound_water) + " [m**3/m**3]" + ", water_init = " + str(self.water_init) + " [m**3/m**3]" + ", Salinity = " + str(self.salinity) + " [mol/L]" + ", Sand content = " + str(self.sand)+ " [m**3/m**3]" + ", Silt content = " + str(self.silt) + " [m**3/m**3]" + ", Clay content = " + str(self.clay) + " [m**3/m**3]" + ", Bulk density = " + str(self.bulk_density) + " [kg/m**3]" + ", Particle density" + str(self.particle_density) + " [kg/m**3]" + ", CEC = " + str(self.CEC) + " [meq/100g]" + ", LOI = " + str(self.LOI) + " [g]" + ", OC = " + str(self.OC) + " [g]" + ", Organic matter = " + str(self.orgm) + " [m**3/m**3]" + ", Bulk real permittivity = " + str(self.bulk_perm) + " [-]" + ", Bulk i permittivity = " + str(self.bulk_iperm) + " [-]" + ", Bulk apparent permittivity  = " + str(self.bulk_aperm) + " [-]" + ', bulk_perm_inf =' + str(self.bulk_perm_inf) + " [-]" + ', bulk_perm_init =' + str(self.bulk_perm_init) + " [-]" + ", Water phase permittivity = " + str(self.water_perm)  + " [-]" + ", Bound water permittivity = " + str(self.bound_water_perm) + " [-]" + ", Water phase i permittivity = " + str(self.water_iperm) + " [-]" + ", Bound water i permittivity = " + str(self.bound_water_iperm) + " [-]" + ", Solid permittivity phase = " + str(self.solid_perm) + " [-]" + ", Air permittivity phase = " + str(self.air_perm) + " [-]" + ", Offset (Hilhorst formula) = " + str(self.offset) + " [-]" + ", Bulk real electr. conductivity = " + str(self.bulk_cond) + " [S/m]" + ", Bulk i electr. conductivity = " + str(self.bulk_icond) + " [S/m]" + ", Water electr. cond = " + str(self.water_cond) + " [S/m]" + ", Matrix electr. conductivity = " + str(self.matrix_cond) + " [S/m]" + "Air electr. conductivity = " + str(self.air_cond) + " [S/m]" + ", Depolarization factor of particles = " + str(self.L) + " [-]" + ", Depolarization factor of water = " + str(self.Lw) + " [-]" + ", Cementation factor (Archie law) = " + str(self.m) + "[-]" + ", Saturation factor (Archie law) = " + str(self.n) + " [-]" + ", Alpha exponent (Roth model) = " + str(self.alpha) + " [-]" + ", Frequency of permittivity measurement = " + str(self.frequency_perm) + " [Hz]" + ", Frequency of electr. conductivity measurement = " + str(self.frequency_cond) + " [Hz]" + " [m]" + ', Land = ' + str(self.land)+ ", Region=" + str(self.region) + ", Depth=" + str(self.depth) + ", Texture=" + str(self.texture) + ', pH_water=' + str(self.pH_water) + ">>>"