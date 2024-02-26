import numpy as np

from .water_from_ec import WaterFromEC
from .water_from_perm import WaterFromPerm
from .frequency_perm import FrequencyPerm
from .frequency_ec import FrequencyEC
from .temperature import Temperature
from .bulk_ec_dc_tc import shift_to_bulk_ec_dc_tc

def Water(soil):
    """
    Return and compute missing values of the soil.df.water attribute using soil.df.bulk_perm or soil.df.bulk_ec.

    This function performs calculations to fill in missing water content values in the soil DataFrame based on bulk permeability (bulk_perm) and electrical conductivity (bulk_ec_dc_tc). 
    It employs different strategies based on available data and sets conditions for water content calculation.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state. Includes columns: water, bulk_perm, frequency_perm, and bulk_ec_dc_tc
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        soil.df.water.values: an array of updated or original soil volumetric water content values

    External functions
    --------
    Temperature : Set missing values of soil.df.temperature and return 
    FrequencyPerm : Set missing values of soil.df.frequency_perm and return 
    WaterFromPerm : Calculate missing values of soil.df.water based on soil.df.bulk_perm
    WaterFromEC : Calculate missing values of soil.df.water based on soil.df.bulk_ec_dc_tc  
    FrequencyEC : Set missing values of soil.df.frequency_ec and return 
    shift_to_bulk_ec_dc_tc : Compute missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec or soil.df.bulk_ec_dc

    Example
    -------
    >>> sample = Soil( bulk_ec = [0.01, 0.02, 0.025, 0.030, 0.040],
                    clay = 10,
                    porosity = 0.47,
                    water_ec = 0.5)

    >>> Water(sample) 
    array([0.105, 0.162, 0.185, 0.206, 0.243])
    """

    # Condition to obtain water from bulk_perm 
    Temperature(soil)
    FrequencyPerm(soil)

    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) for x in range(soil.n_states)) and (np.isnan(soil.df.frequency_perm)).all():
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Unmodified value. Please provide soil.frequency_perm" if True else soil.info.bulk_perm[x] for x in range(soil.n_states)]
    
    elif any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) for x in range(soil.n_states)) and not (np.isnan(soil.df.frequency_perm).all()):
        WaterFromPerm(soil) 

    # Condition to obtain water from bulk_ec_dc_tc
    FrequencyEC(soil)
    shift_to_bulk_ec_dc_tc(soil)
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
        WaterFromEC(soil)        

    # Converting negative results due to fitting to zero
    soil.info['water'] = [str(soil.info.water[x]) + "--> Set to 0 because of < 0 results" if soil.df.water[x]<0 
                          or soil.info.water[x] ==str(soil.info.water[x]) + "--> Set to 0 because of < 0 results"
                            else soil.info.water[x] for x in range(soil.n_states)]

    soil.df['water'] = [ 0 if soil.df.water[x]<0 else soil.df.water[x] for x in range(soil.n_states)] 

    return soil.df.water.values