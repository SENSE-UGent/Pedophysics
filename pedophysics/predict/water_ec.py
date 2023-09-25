import numpy as np
from scipy.optimize import minimize, Bounds

from pedophysics.pedophysical_models.water_ec import SenGoode
from pedophysics.pedophysical_models.bulk_ec import Fu, Rhoades
from pedophysics.pedophysical_models.bulk_perm import Hilhorst
from pedophysics.utils.stats import R2_score

from .temperature import Temperature
from .particle_density import ParticleDensity
from .solid_ec import SolidEC
from .texture import Texture
from .frequency_ec import FrequencyEC
from .water_perm import WaterPerm


def WaterEC(soil):
    """
    Return and computes soil.df.water_ec based on soil.df.water, soil.df.bulk_ec and soil.df.bulk_perm. 

    This function uses either a non-fitting or a fitting approach based on the available soil 
    attributes and conditions to estimate the water EC. The non-fitting approach is employed 
    when either the water EC values are missing and soil salinity values are available, or only 
    one state has non-missing values for water and bulk EC. Otherwise, if there are at least two 
    states with non-missing values for bulk EC and either water or bulk permittivity, the fitting 
    approach is used.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water_ec, salinity, water, bulk_ec, frequency_ec, and bulk_perm. 
        - n_states : int
            Number of soil states.

    Returns
    -------
    numpy.ndarray
        Array containing the updated soil water pore real electrical conductivity values. 

    Notes
    -----
    The function first ensures that the temperature and frequency EC attributes are available.
    Then, based on the conditions specified, it employs the suitable approach to estimate the water EC.

    External Functions
    ------------------
    - Temperature : Function to ensure the temperature values are available for the soil.
    - FrequencyEC : Function to ensure the frequency EC values are available for the soil.
    - non_fitting : Non-fitting approach function for estimating water EC.
    - fitting     : Fitting approach function for estimating water EC.

    Example
    -------
    >>> sample = Soil( bulk_ec=[0.02, 0.03, 0.04, 0.05, 0.06], 
                bulk_perm=[11.5, 14.8, 17, 20, 22.7],
                clay=5,
                bulk_density=1.48,
                instrument='TDR')

    >>> predict.WaterEC(sample)
    array([0.289855, 0.289855, 0.289855, 0.289855, 0.289855])
    """

    Temperature(soil)
    FrequencyEC(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.salinity[x]) for x in range(soil.n_states)) or sum(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)) == 1:
        non_fitting(soil)

    # Condition for fitting approach
    if sum(not np.isnan(soil.bulk_ec[x]) and (not np.isnan(soil.water[x]) or not np.isnan(soil.bulk_perm[x])) and np.isnan(soil.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting(soil)

    return soil.df.water_ec.values


def non_fitting(soil):
    """
    Decides between two different non-fitting approaches to obtain soil.df.water_ec. 

    This function employs two non-fitting strategies based on the available data in the soil object:
    1. Estimates the water EC using the salinity attribute if water EC values are missing but salinity values are available.
    2. Estimates the water EC from the bulk EC attribute if water EC values are missing but values for water and bulk EC are available.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water_ec, salinity, water, bulk_ec, frequency_ec, and bulk_perm.
        - n_states : int
            Number of soil states.

    External Functions
    ------------------
    - from_salinity : Function that estimates water EC from the salinity attribute.
    - from_bulk_ec  : Function that estimates water EC from the bulk EC attribute.
    """

    # Condition for non-fitting approach using salinity
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.salinity[x]) for x in range(soil.n_states)):
        from_salinity(soil)

    # Condition for non-fitting approach using bulk_ec
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)):
        from_bulk_ec(soil)


def from_salinity(soil):
    """
    Calculates soil.df.water_ec based on soil.df.salinity using the SenGoode function.

    This function uses the salinity attribute in the soil object to estimate missing water EC values. The 
    SenGoode function is used to perform this conversion. For each state in the soil object, if the water EC 
    value is missing but the salinity value is available, the water EC is calculated using the SenGoode function.
    
    Additionally, an information note is added to the `info` attribute of the soil object indicating the method 
    of calculation for the water EC.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water_ec, salinity, and temperature.
        - n_states : int
            Number of soil states.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    External Functions
    ------------------
    - SenGoode : Function that estimates water EC from the temperature and salinity attributes.
    """

    # Calculating and saving water_ec and its info
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated using SenGood function in predict.water_ec.from_salinity" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated using SenGood function in predict.water_ec.from_salinity"
                                else soil.info.water_ec[x] for x in range(soil.n_states)]
 
    soil.df['water_ec'] = [SenGoode(soil.df.temperature[x], soil.df.salinity[x]) if np.isnan(soil.df.water_ec[x]) 
                                else soil.df.water_ec[x] for x in range(soil.n_states)]


def from_bulk_ec(soil):
    """
    Calculates soil.df.water_ec based on soil.df.bulk_ec using the Fu function.

    This function uses the bulk EC and other soil properties to estimate missing water EC values. The 
    Fu function, which is a pedophysical model for bulk electrical conductivity, is used to perform the estimation.
    For each state in the soil object, if the water EC value is missing, a minimization routine is applied 
    to determine the optimal water EC that best fits the observed bulk EC.

    Additionally, an information note is added to the `info` attribute of the soil object indicating the method 
    of calculation for the water EC using the Fu function.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water_ec, clay, bulk_density, particle_density, solid_ec, dry_ec, sat_ec, and bulk_ec.
        - n_states : int
            Number of soil states.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    External Functions
    ------------------
    - Texture : Function that provides the texture attributes of the soil.
    - ParticleDensity : Function that provides the particle density of the soil.
    - SolidEC : Function that provides the electrical conductivity of the solid particles.
    - Fu : Model function to estimate bulk EC from various soil properties.
    """

    Texture(soil)
    ParticleDensity(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water_ec
    def objective_wat_ec(water_ec, wat, clay, bulk_density, particle_density, solid_ec, dry_ec, sat_ec, EC):
        return abs(Fu(wat, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec) - EC)
    
    # Calculating optimal water_ec
    wat_ec = []
    for i in range(soil.n_states):
        res = minimize(objective_wat_ec, 0.14, args=(soil.df.water[i], soil.df.clay[i], soil.df.bulk_density[i], soil.df.particle_density[i], soil.df.solid_ec[i], 
                                                     soil.df.dry_ec[i], soil.df.sat_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 2)] )
        wat_ec.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

    # Saving calculated water_ec and its info
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_ec.from_bulk_ec" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_ec.from_bulk_ec"
                                 else soil.info.water_ec[x] for x in range(soil.n_states)]

    soil.df['water_ec'] = [round(wat_ec[x], soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water[x] for x in range(soil.n_states) ]


def fitting(soil):
    """
    Decides between two different fitting approaches to obtain soil.df.water_ec. 

    Depending on the available data within the soil object, this function selects one of two fitting 
    functions: Rhoades function or Hilhorst function. The selection criterion is based on the number of 
    missing water EC values and the availability of certain soil properties (e.g., water content and bulk permittivity).
    If enough data points are available and the conditions are met, the fitting function is applied to estimate the 
    missing water EC values.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_ec, water, bulk_perm, and water_ec.
        - n_states : int
            Number of soil states.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    External Functions
    ------------------
    - fitting_rhoades : Function that applies the Rhoades function to estimate missing water EC values using bulk EC and water content.
    - fitting_hilhorst : Function that applies the Hilhorst function to estimate missing water EC values using bulk EC and bulk permittivity.
    """

    # Condition for fitting approach using Rhoades function
    if sum(not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting_rhoades(soil)
    
    # Condition for fitting approach using Rhoades function
    elif sum(not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting_hilhorst(soil)


def fitting_rhoades(soil):
    """
    Calculates soil.df.water_ec using the Rhoades function based on soil properties.

    The function estimates water's EC by fitting the Rhoades function to the available bulk EC and water content data.
    In addition to estimating the missing water EC values, this function calculates other parameters such as s_ec, E, 
    and F that are essential for the Rhoades function.
    The function initially selects calibration data by identifying valid data points from the soil.df attributes. 
    Then, it utilizes the `minimize` method to fit the Rhoades function, allowing for the estimation of the required 
    parameters. Finally, the calculated parameters are stored in the soil object, and the fitting's R2 score is recorded 
    in the info attribute of the soil object.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_ec, water, s_ec, and water_ec.
        - n_states : int
            Number of soil states.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - E : float
            An attribute in the soil object where the calculated E value is stored.
        - F : float
            An attribute in the soil object where the calculated F value is stored.

    External Functions
    ------------------
    - Rhoades : A mathematical function used to relate bulk EC with water content, water EC, and other parameters.
    """

    # Selecting calibration data
    arg_EC_wn = np.array([soil.df.bulk_ec[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    arg_water_wn = np.array([soil.df.water[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    
    # Removing NaNs from calibration data
    valid_indices = ~np.isnan(arg_EC_wn) & ~np.isnan(arg_water_wn)
    arg_EC = arg_EC_wn[valid_indices]
    arg_water = arg_water_wn[valid_indices]
    
    # Define the initial guesses
    bounds = Bounds([0.00001, 0], [2, 0.1])
    initial_guess_watec = 0.15
    initial_guess_s_ec = 0
    initial_guess_E = 1
    initial_guess_F = 0.38

    # Defining minimization function to obtain water_ec and s_ec while fixing E and F
    def objective_water_ec(params, wat, bulk_ec, E, F):
        water_ec, s_ec = params
        residuals = (Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec)**2
        return np.sum(residuals)

    # Calculating optimal water_ec and s_ec
    res1 = minimize(objective_water_ec, [initial_guess_watec, initial_guess_s_ec], args=(arg_water, arg_EC, initial_guess_E, initial_guess_F), bounds=bounds)
    best_water_ec, best_s_ecs = res1.x

    # Saving calculated s_ec and its info
    soil.info['s_ec'] = [str(soil.info.s_ec[x]) + "--> Calculated by fitting Rhoades function in predict.water_ec.fitting_rhoades" if np.isnan(soil.df.s_ec[x])
                            or soil.info.s_ec[x] == str(soil.info.s_ec[x]) + "--> Calculated by fitting Rhoades function in predict.water_ec.fitting_rhoades"
                            else soil.info.s_ec[x] for x in range(soil.n_states)]
    
    soil.df['s_ec'] = [round(best_s_ecs, soil.roundn+3) if np.isnan(soil.df.s_ec[x]) else soil.df.s_ec[x] for x in range(soil.n_states) ]

    # Defining minimization function to obtain E and F while fixing water_ec and s_ec
    def objective_others(params, wat, bulk_ec, water_ec, s_ec):
        E, F = params
        residuals = np.sum((Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec)**2)
        return residuals

    # Calculating optimal E and F
    res2 = minimize(objective_others, [initial_guess_E, initial_guess_F], args=(arg_water, arg_EC, best_water_ec, best_s_ecs))
    best_E, best_F = res2.x
    soil.E = best_E
    soil.F = best_F

    # Calculating the R2 score of the fitting
    R2 = round(R2_score(arg_EC, Rhoades(arg_water, best_water_ec, best_s_ecs, best_E, best_F)), soil.roundn)
    
    # Saving calculated water_ec and its info with R2
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2 = "+str(R2)+") Rhoades function in predict.water_ec.fitting_rhoades" if np.isnan(soil.df.water_ec[x]) 
                            or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2 = "+str(R2)+") Rhoades function in predict.water_ec.fitting_rhoades"
                            else soil.info.water_ec[x] for x in range(soil.n_states)]

    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water_ec[x] for x in range(soil.n_states) ]


def fitting_hilhorst(soil):
    """
    Calculates soil.df.water_ec using the Hilhorst function based on soil properties.
    
    This function estimates water's EC by fitting the Hilhorst function to the available bulk EC, 
    bulk permittivity and water permittivity data. Additionally, the function calculates the 
    `offset_perm` parameter essential for the Hilhorst function.
    
    The function begins by selecting calibration data by identifying valid data points from the soil.df 
    attributes, and subsequently applies a filtering mechanism to eliminate NaN values. The `minimize` 
    method is then used to fit the Hilhorst function, enabling the estimation of the desired parameters. 
    Finally, the derived parameters are stored within the soil object, and the fitting's R2 score is 
    added to the info attribute of the soil object.
    
    Parameters
    ----------
    soil : object
        A custom soil object containing:
        
        df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            This includes attributes such as bulk_ec, bulk_perm, water_perm, offset_perm, and water_ec.
        n_states : int
            The count of soil states.
        info : dict
            A dictionary which conveys descriptive data on how each array-like attribute was either determined or modified.
        roundn : int
            Number of decimal places to round results.

    External Functions
    ------------------
    WaterPerm : function
        Used for deriving the water permittivity for each soil state.
    Hilhorst : function
        A mathematical equation to relate bulk EC with water content, water EC, water permittivity, and other parameters.
    """
    
    WaterPerm(soil)

    # Selecting calibration data
    arg_EC_wn = np.array([soil.df.bulk_ec[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                            else np.nan for x in range(soil.n_states)])
    arg_bulk_perm_wn = np.array([soil.df.bulk_perm[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                              else np.nan for x in range(soil.n_states)])
    arg_water_perm_wn = np.array([soil.df.water_perm[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                               else np.nan for x in range(soil.n_states)])

    # Removing NaNs from calibration data
    valid_indices = ~np.isnan(arg_EC_wn) & ~np.isnan(arg_bulk_perm_wn)
    arg_EC = arg_EC_wn[valid_indices]
    arg_bulk_perm = arg_bulk_perm_wn[valid_indices]
    arg_water_perm = arg_water_perm_wn[valid_indices]
    
    # Define the initial guesses
    bounds = Bounds([0.00001, -10], [2, 10])
    initial_guess_offset_perm = 4
    initial_guess_watec = 0.15

    # Defining minimization function
    def objective_water_ec(param, bulk_perm, bulk_ec, water_perm):
        water_ec, offset_perm = param
        residuals = (Hilhorst(bulk_ec, water_ec, water_perm, offset_perm) - bulk_perm)**2
        return np.sum(residuals)

    # Calculating optimal water_ec and offset_perm
    res = minimize(objective_water_ec, [initial_guess_watec, initial_guess_offset_perm], args=(arg_bulk_perm, arg_EC, arg_water_perm), bounds=bounds)
    best_water_ec, best_offset_perm = res.x

    # Saving calculated offset_perm and its info
    soil.info['offset_perm'] = [str(soil.info.offset_perm[x]) + "--> Calculated by fitting Hilhorst function in predict.water_ec.fitting_hilhorst" if np.isnan(soil.df.offset_perm[x]) 
                                or soil.info.offset_perm[x] == str(soil.info.offset_perm[x]) + "--> Calculated by fitting Hilhorst function in predict.water_ec.fitting_hilhorst"
                                 else soil.info.offset_perm[x] for x in range(soil.n_states)]
    
    soil.df['offset_perm'] = [round(best_offset_perm, soil.roundn+3) if np.isnan(soil.df.offset_perm[x]) else soil.df.offset_perm[x] for x in range(soil.n_states) ]

    # Calculating the R2 score of the fitting
    R2 = round(R2_score(arg_bulk_perm, Hilhorst(arg_EC, best_water_ec, arg_water_perm, best_offset_perm)), soil.roundn)
    
    # Saving calculated water_ec and its info with R2
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") Hilhorst function in predict.water_ec.fitting_hilhorst" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") Hilhorst function in predict.water_ec.fitting_hilhorst"
                                 else soil.info.water_ec[x] for x in range(soil.n_states)]
    
    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water_ec[x] for x in range(soil.n_states) ]