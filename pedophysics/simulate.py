import numpy as np
import pandas as pd

class Soil(object):
    """
    A class to represent a soil sample with its characteristics.

    Attributes
    ----------
    temperature : array-like
        Soil bulk temperature [K]
    water : array-like
        Soil volumetric water content [m**3/m**3]
    salinity : array-like
        Soil salinity (NaCl) of the bulk pore fluid [mol/L]
    sand : array-like
        Soil sand content [g/g]*100
    silt : array-like
        Soil silt content [g/g]*100
    clay : array-like
        Soil clay content [g/g]*100
    porosity : array-like
        Soil porosity [m**3/m**3]
    bulk_density : array-like 
        Soil bulk density [kg/m**3]
    particle_density : array-like
        Soil particle density [kg/m**3]
    CEC : array-like
        Soil cation exchange capacity [meq/100g]
    orgm : array-like
        Soil organic matter [g/g]*100
    bulk_perm : array-like
        Soil bulk real relative dielectric permittivity [-]
    bulk_perm_inf : array-like
        Soil bulk real relative permittivity at infinite frequency [-]
    water_perm : array-like
        Soil water phase real dielectric permittivity [-]
    solid_perm : array-like
        Soil solid real relative dielectric permittivity phase [-]
    air_perm : array-like
        Soil air real relative dielectric permittivity phase [-]
    offset_perm : array-like
        Soil bulk real relative dielectric permittivity when soil bulk real electrical conductivity is zero [-]
    bulk_ec : array-like
        Soil bulk real electrical conductivity [S/m]
    bulk_ec_dc : array-like
        Soil bulk real electrical conductivity direct current [S/m]
    bulk_ec_dc_tc : array-like
        Soil bulk real electrical conductivity direct current (0 Hz) temperature corrected (298.15 K) [S/m]
    water_ec : array-like
        Soil water real electrical conductivity [S/m]
    s_ec : array-like
        Soil bulk real surface electrical conductivity [S/m]
    solid_ec : array-like
        Soil solid real electrical conductivity [S/m]
    dry_ec : array-like
        Soil bulk real electrical conductivity at zero water content [S/m]
    sat_ec : array-like 
        Soil bulk real electrical conductivity at saturation water content [S/m]
    frequency_perm : array-like
        Frequency of dielectric permittivity measurement [Hz]
    frequency_ec : array-like
        Frequency of electric conductivity measurement [Hz]
    L : single-value
        Soil scalar depolarization factor of solid particles (effective medium theory) [-]
    Lw : single-value 
        Soil scalar depolarization factor of water aggregates (effective medium theory) [-]
    m : single-value
        Soil cementation factor as defined in Archie law [-]
    n : single-value
        Soil saturation factor as defined in Archie second law [-]
    alpha : single-value
        Soil alpha exponent as defined in volumetric mixing theory [-]
    texture : str
        Soil texture according to USDA convention: "Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Silty clay loam", "Sandy clay", "Clay", "Silty clay"
    instrument : str
        Instrument utilized: 'HydraProbe', 'TDR', 'GPR', 'Miller 400D', 'Dualem'
    info : DataFrame
        Data Frame containing descriptive information about how each array-like attribute was determined or modified.
    df : DataFrame
        Data Frame containing the quantitative information of all soil array-like attributes for each state.
    E : single-value
        Empirical constant as in Rohades model [-]
    F : single-value
        Empirical constant as in Rohades model [-]
    roundn : int
        Number of decimal places to round results.
    range_ratio : single-value
        Factor for extending extrapolation domain during fitting modelling
    n_states : int
        Number of soil states

    Notes
    -----
    Attributes provided by the user that do not match the expected types or values 
    will raise a ValueError.
    """

    def __init__(self, **kwargs):
        # Define acceptable types for each argument
        array_like_types = [float, np.float64, int, list, np.ndarray]
        single_value = [float, np.float64, int]
        attributes = {
                'temperature': array_like_types,
                'water': array_like_types,
                'salinity': array_like_types,
                'sand': array_like_types,
                'silt': array_like_types,
                'clay': array_like_types,
                'porosity': array_like_types,
                'bulk_density': array_like_types,
                'particle_density': array_like_types,
                'CEC': array_like_types,
                'orgm': array_like_types,
                'bulk_perm': array_like_types,
                'bulk_perm_inf': array_like_types,
                'air_perm': array_like_types,
                'water_perm': array_like_types,
                'solid_perm': array_like_types,
                'offset_perm': array_like_types,
                'bulk_ec': array_like_types,
                'bulk_ec_tc': array_like_types,
                'bulk_ec_dc': array_like_types,
                'bulk_ec_dc_tc': array_like_types,
                'water_ec': array_like_types,
                'solid_ec': array_like_types,
                'dry_ec': array_like_types,
                'sat_ec': array_like_types,
                's_ec': array_like_types,
                'frequency_perm': array_like_types,
                'frequency_ec': array_like_types,
                'L': single_value,
                'Lw': single_value,
                'm': single_value,
                'n': single_value,
                'alpha': single_value,
                'texture': [str],
                'instrument': [str],
                'range_ratio': single_value,
                'n_states': single_value,
                'E': single_value,
                'F': single_value,
                'roundn': [int]
                }

        accepted_values = {
            'texture': ["Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Sandy clay", "Clay", "Silty clay", np.nan],
            'instrument': ["TDR", "GPR", 'HydraProbe', 'EMI Dualem', 'EMI EM38-DD', np.nan]
        }

        # Convert all inputs to np.ndarray if they are of type list, int, or float
        def to_ndarray(arg, key=None):
            if key in ['texture', 'instrument']:
                return arg  # return the argument if it is 'texture' or 'instrument'
            if isinstance(arg, (list, int, np.float64, float)):
                return np.array([arg]) if isinstance(arg, (int, np.float64, float)) else np.array(arg)
            return arg

        # Check each input argument
        for key in attributes:

            if key in kwargs:
                value = kwargs[key]

                if type(value) in attributes[key]:
                    # if the key is 'texture' or 'instrument' verify if value is in the accepted_values
                    if key in ['texture', 'instrument'] and value not in accepted_values[key]:
                        raise ValueError(f"Invalid value for '{key}'. Must be one of {accepted_values[key]}")
                    setattr(self, key, to_ndarray(value, key=key))
                else:
                    raise ValueError(f"'{key}' must be one of {attributes[key]}")
                
            else:
                # If the key is not provided in the kwargs, set it as np.nan.
                setattr(self, key, to_ndarray(np.nan, key=key))
            
        self.roundn = 3 if np.isnan(self.roundn[0]) else self.roundn
        self.range_ratio = 2 if np.isnan(self.range_ratio[0]) else self.range_ratio
        
        ### Fill the state variables with nans when are shorter than n_states
        array_like_attributes = ['temperature', 'water', 'salinity', 'sand', 'silt', 'clay', 'porosity', 'bulk_density', 'particle_density', 'CEC',
                            'orgm', 'bulk_perm', 'bulk_perm_inf', 'air_perm', 'water_perm', 'solid_perm', 'offset_perm', 
                            'bulk_ec', 'bulk_ec_tc', 'bulk_ec_dc', 'bulk_ec_dc_tc', 'water_ec', 'solid_ec', 'dry_ec', 'sat_ec', 's_ec', 'frequency_perm', 'frequency_ec']

        # calculate the max length of the input arrays
        n_states = max([len(getattr(self, attr)) for attr in array_like_attributes])
        self.n_states = n_states                            # Number of states of the soil

        # Now loop over each attribute in the list
        for attribute in array_like_attributes:
            attr = getattr(self, attribute)
            
            if len(attr) != n_states:
                setattr(self, attribute, np.append(attr, [np.nan]*(n_states - len(attr))))        

            if ~np.isnan(attr[0]) and (np.isnan(attr[1:(n_states)])).all():
                setattr(self, attribute, np.append(attr[0], [attr[0]]*(n_states - 1)))     

        ### Defining special attributes ### 
        self.df = pd.DataFrame({attr: getattr(self, attr) for attr in array_like_attributes})

        # defining soil.info
        self.info = self.df.where(pd.notna(self.df), 'NaN')
        self.info = self.info.where(pd.isna(self.info), 'Value given by the user')
        
    # Simplify the getter methods using __getattr__
    def __getattr__(self, name):
        """
        Custom attribute access mechanism.

        Parameters
        ----------
        name : str
            Name of the attribute to be accessed.

        Returns
        -------
        np.ndarray
            The value of the attribute.

        Raises
        ------
        AttributeError
            If the attribute does not exist.
        """
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError(f"No such attribute: {name}")        
        
    def __str__(self):                                                   
        """
        Return a string representation of the class.

        Returns
        -------
        str
            String representation of the class as Soil.df
        """
        return str(self.df)