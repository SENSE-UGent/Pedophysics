import numpy as np
import pandas as pd

class Soil(object):
    def __init__(self, **kwargs):
        # Define acceptable types for each argument

        attributes = {
                'temperature': [float, np.ndarray, int, list],
                'water': [float, np.ndarray, int, list],
                'salinity': [float, np.ndarray, int, list],
                'sand': [float, int, np.ndarray, list],
                'silt': [float, int, np.ndarray, list],
                'clay': [float, int, np.ndarray, list],
                'bulk_density': [float, np.ndarray, int, list],
                'particle_density': [float, np.ndarray, int, list],
                'CEC': [float, np.ndarray, int, list],
                'OC': [float, np.ndarray, int, list],
                'orgm': [float, np.ndarray, int, list],
                'bulk_perm': [float, np.ndarray, int, list],
                'bulk_perm_inf': [float, np.ndarray, int, list],
                'air_perm': [float, np.ndarray, int, list],
                'water_perm': [float, np.ndarray, int, list],
                'solid_perm': [float, np.ndarray, int, list],
                'offset_perm': [float, np.ndarray, int, list],
                'bulk_ec': [float, np.ndarray, int, list],
                'water_ec': [float, np.ndarray, int, list],
                'solid_ec': [float, np.ndarray, int, list],
                'dry_ec': [float, np.ndarray, int, list],
                'sat_ec': [float, np.ndarray, int, list],
                's_ec': [float, np.ndarray, int, list],
                'frequency_perm': [float, np.ndarray, int, list],
                'frequency_ec': [float, np.ndarray, int, list],
                'L': [float, int],
                'Lw': [float, int],
                'm': [float, int],
                'n': [float, int],
                'alpha': [float, int],
                'E': [float, int],
                'F': [float, int],
                'range_ratio': [float, int],
                'n_states': [float, int],
                'texture': [str],
                'instrument': [str],
                'roundn': [int],
                }

        accepted_values = {
            'texture': ["Sand", "Loamy sand", "Sandy loam", "Loam", "Silt loam", "Silt", "Sandy clay loam", "Clay loam", "Sandy clay", "Clay", "Silty clay", np.nan],
            'instrument': ["TDR", "GPR", 'HydraProbe', 'EMI Dualem', 'EMI EM38-DD', np.nan]
        }

        # Convert all inputs to np.ndarray if they are of type list, int, or float
        def to_ndarray(arg, key=None):
            if key in ['texture', 'instrument']:
                return arg  # return the argument if it is 'texture' or 'instrument'
            if isinstance(arg, (list, int, float)):
                return np.array([arg]) if isinstance(arg, (int, float)) else np.array(arg)
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
        state_attribute = ['temperature', 'water', 'salinity', 'sand', 'silt', 'clay', 'bulk_density', 'particle_density', 'CEC', 'OC',
                'orgm', 'bulk_perm', 'bulk_perm_inf', 'air_perm', 'water_perm', 'solid_perm', 'offset_perm', 'bulk_ec', 'water_ec',
                'solid_ec', 'dry_ec', 'sat_ec', 's_ec', 'frequency_perm', 'frequency_ec']

        # calculate the max length of the input arrays
        n_states = max([len(getattr(self, attr)) for attr in state_attribute])
        self.n_states = n_states                            # Number of states of the soil

        # Now loop over each attribute in the list
        for attribute in state_attribute:
            attr = getattr(self, attribute)
            
            if len(attr) != n_states:
                setattr(self, attribute, np.append(attr, [np.nan]*(n_states - len(attr))))        

            if ~np.isnan(attr[0]) and (np.isnan(attr[1:(n_states)])).all():
                setattr(self, attribute, np.append(attr[0], [attr[0]]*(n_states - 1)))     

        ### Defining Soil.df and Soil.info ##### 
        self.df = pd.DataFrame({attr: getattr(self, attr) for attr in state_attribute})

        # defining soil.info
        self.info = self.df.where(pd.notna(self.df), np.nan)
        self.info = self.info.where(pd.isna(self.info), 'Values given by the user')
        
    # Simplify the getter methods using __getattr__
    def __getattr__(self, name):
        if name in self.__dict__:
            return self.__dict__[name]
        else:
            raise AttributeError(f"No such attribute: {name}")        
        
    def __str__(self):                                                   
        """ Returns a string representation of self """
        return str(self.df)