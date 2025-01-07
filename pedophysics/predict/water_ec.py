import numpy as np
from scipy.optimize import minimize, Bounds

from pedophysics.pedophysical_models.water_ec import SenGoode
from pedophysics.pedophysical_models.bulk_ec import Fu, Rhoades
from pedophysics.pedophysical_models.bulk_perm import Hilhorst
from pedophysics.utils.stats import R2_score

from .temperature import Temperature
from .porosity import Porosity
from .solid_ec import SolidEC
from .texture import Texture
from .frequency_ec import FrequencyEC
from .water_perm import WaterPerm

from .bulk_ec_dc_tc import shift_to_bulk_ec_dc_tc

def WaterEC(soil):
    """
    Compute missing values of soil.df.water_ec and return  

    This function evaluates soil states to determine the appropriate approach for estimating water EC. 
    It considers non-fitting approaches based on salinity and bulk electrical conductivity, 
    as well as fitting approaches using the Rhoades or Hilhorst functions depending on the availability of water content, bulk electrical conductivity, and bulk permeability data.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: `water`, `bulk_ec_dc_tc`, `water_ec`, `salinity`, and potentially `bulk_perm`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - n_states: int
            The number of soil states represented in the `df`.

    Returns
    -------
    numpy.ndarray
        soil.df.water_ec.values: array containing the updated soil water pore real electrical conductivity values. 

    Notes
    -----
    - Non-fitting approaches are applied when specific data are available, such as salinity or bulk electrical conductivity, without the need for additional parameters.
    - Fitting approaches, such as those using the Rhoades or Hilhorst functions, are applied when there are sufficient data points with known water content or bulk permeability.

    External Functions
    ------------------
    - Temperature : Set missing values of soil.df.temperature and return 
    - FrequencyEC : Set missing values of soil.df.frequency_ec and return 
    - shift_to_bulk_ec_dc_tc : Compute missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec or soil.df.bulk_ec_dc
    - from_salinity	: Calculate missing values of soil.df.water_ec based on soil.df.salinity 
    - from_ec : Calculate missing values of soil.df.water_ec based on soil.df.bulk_ec_dc_tc
    - fitting_rhoades : Calculate missing values of soil.df.water_ec using the Rhoades function in a fitting approach
    - fitting_hiolhorst : Calculate missing values of soil.df.water_ec using the Hilhorst function in a fitting approach

    Example
    -------
    >>> sample = Soil( bulk_ec=[0.02, 0.03, 0.04, 0.05, 0.06], 
                bulk_perm=[11.5, 14.8, 17, 20, 22.7],
                clay=5,
                porosity=0.44,
                instrument='TDR')

    >>> predict.WaterEC(sample)
    array([0.289855, 0.289855, 0.289855, 0.289855, 0.289855])
    """
    Temperature(soil)
    FrequencyEC(soil)
    shift_to_bulk_ec_dc_tc(soil)

    # Condition for non-fitting approach using salinity
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.salinity[x]) for x in range(soil.n_states)):
        from_salinity(soil)

    # Conditions for fitting approaches
    if sum(not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.water[x]) and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2 or sum(not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:

        # Condition for fitting approach using Rhoades function
        if sum(not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.water[x]) and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
            fitting_rhoades(soil)
        
        # Condition for fitting approach using Rhoades function
        elif sum(not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
            fitting_hilhorst(soil)

    # Condition for non-fitting approach using bulk_ec_dc_tc
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
        from_ec(soil)

    return soil.df.water_ec.values


def from_salinity(soil):
    """
    Calculate missing values of soil.df.water_ec based on soil.df.salinity 

    This function estimates water EC for each soil state based on temperature and salinity data using the SenGood function. 
    Estimated water EC values are saved in the soil DataFrame, and an information string indicating the use of the SenGood function 
    for calculation is appended to each relevant entry in the soil information dictionary.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: `temperature`, `salinity`, and potentially `water_ec`.
        - n_states: int
            The number of soil states represented in the `df`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes with the estimated water EC values and does not return any value.

    External Functions
    ------------------
    - SenGoode : Calculate soil water real electrical conductivity using the Sen and Goode model and return
    """

    # Check for missing values
    missing_water_ec_before = soil.df['water_ec'].isna()

    soil.df['water_ec'] = [SenGoode(soil.df.temperature[x], soil.df.salinity[x]) 
                            if np.isnan(soil.df.water_ec[x]) 
                            else soil.df.water_ec[x] for x in range(soil.n_states)]
    missing_water_ec_after = soil.df['water_ec'].isna()
    
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + (
            "--> Calculated using SenGood function in predict.water_ec.from_salinity"
            if missing_water_ec_before[x] and not missing_water_ec_after[x]
            #else "--> Provide water_ec, otherwise salinity"
            #if missing_water_ec_before[x] and missing_water_ec_after[x]
            else "")
        if missing_water_ec_before[x]
        else soil.info.water_ec[x]
        for x in range(soil.n_states)]


def from_ec(soil):
    """
    Calculate missing values of soil.df.water_ec based on soil.df.bulk_ec_dc_tc

    This function applies the Fu function within a minimization process to estimate soil water EC based on soil properties such as 
    water content, clay content, porosity, solid EC, dry EC, and saturated EC. The estimation is performed for each soil state where water EC is unknown.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: `water`, `clay`, `porosity`, `solid_ec`, `dry_ec`, `sat_ec`, `bulk_ec_dc_tc`, and potentially `water_ec`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - n_states: int
            The number of soil states represented in the `df`.
        - roundn: int
            The number of decimal places for rounding estimated water EC values.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes with the estimated water EC values and does not return any value.

    External Functions
    ------------------
    - Texture : Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    - Porosity : Calculate missing values of soil.df.porosity and return
    - SolidEC : Set missing values of soil.df.solid_ec and return
    - Fu : Calculate the soil bulk real electrical conductivity using the Fu model and return

    """
    Texture(soil)
    Porosity(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water_ec
    def objective_wat_ec(water_ec, wat, clay, porosity, solid_ec, dry_ec, sat_ec, EC):
        return abs(Fu(wat, clay, porosity, water_ec, solid_ec, dry_ec, sat_ec) - EC)
    
    # Calculating optimal water_ec
    wat_ec = []
    for i in range(soil.n_states):
        res = minimize(objective_wat_ec, 0.14, args=(soil.df.water[i], soil.df.clay[i], soil.df.porosity[i], soil.df.solid_ec[i], 
                                                     soil.df.dry_ec[i], soil.df.sat_ec[i], soil.df.bulk_ec_dc_tc[i]), bounds=[(0, 2)] )
        wat_ec.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

    # Saving calculated water_ec and its info
    missing_water_ec_before = soil.df['water_ec'].isna()

    soil.df['water_ec'] = [round(wat_ec[x], soil.roundn+3) 
                           if np.isnan(soil.df.water_ec[x]) 
                           else soil.df.water[x] for x in range(soil.n_states) ]
    
    missing_water_ec_after = soil.df['water_ec'].isna()
    
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + (
            "--> Calculated using Fu function (reported R2=0.98) in predict.water_ec.from_ec"
            if missing_water_ec_before[x] and not missing_water_ec_after[x]
            else "--> Provide water_ec; otherwise bulk_ec_dc_tc, water, clay and porosity"
            if missing_water_ec_before[x] and missing_water_ec_after[x]
            else "")
        if missing_water_ec_before[x]
        else soil.info.water_ec[x]
        for x in range(soil.n_states)]
    

def fitting_rhoades(soil):
    """
    Calculate missing values of soil.df.water_ec using the Rhoades function in a fitting approach

    This function selects calibration data based on available water content and bulk electrical conductivity data, removes NaNs, 
    and uses the Rhoades function within a minimization process to estimate `water_ec` and `s_ec` parameters. 
    It then fixes these parameters to estimate the remaining parameters of the Rhoades function, `E` and `F`. The quality of the fit is evaluated using the R2 score.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            includes: `water`, `bulk_ec_dc_tc`, `water_ec`, `s_ec`, and potentially other related parameters.
        - n_states: int
            The number of soil states represented in the `df`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - roundn: int
            The number of decimal places for rounding estimated parameter values.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes with the estimated parameters and does not return any value.

    Notes
    -----
    - The fitting process involves two steps: first, estimating `water_ec` and `s_ec` with fixed `E` and `F`, and second, estimating `E` and `F` with fixed `water_ec` and `s_ec`.
    - The process uses calibration data where both water content and bulk electrical conductivity are known.

    External Functions
    ------------------
    - Rhoades : Calculate the soil bulk real electrical conductivity using the Rhoades model and return
    - R2_score : Calculate the coefficient of determination (R^2) of a prediction and return.
    """
    # Selecting calibration data

    arg_EC_wn = np.array([soil.df.bulk_ec_dc_tc[x] if not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    arg_water_wn = np.array([soil.df.water[x] if not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    
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
    def objective_water_ec(params, wat, bulk_ec_dc_tc, E, F):
        water_ec, s_ec = params
        residuals = (Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec_dc_tc)**2
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
    def objective_others(params, wat, bulk_ec_dc_tc, water_ec, s_ec):
        E, F = params
        residuals = np.sum((Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec_dc_tc)**2)
        return residuals

    # Calculating optimal E and F
    res2 = minimize(objective_others, [initial_guess_E, initial_guess_F], args=(arg_water, arg_EC, best_water_ec, best_s_ecs))
    best_E, best_F = res2.x
    soil.E = best_E
    soil.F = best_F

    # Calculating the R2 score of the fitting
    R2 = round(R2_score(arg_EC, Rhoades(arg_water, best_water_ec, best_s_ecs, best_E, best_F)), soil.roundn)

    missing_water_ec_before = soil.df['water_ec'].isna()

    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) 
                           if np.isnan(soil.df.water_ec[x]) 
                           else soil.df.water_ec[x] for x in range(soil.n_states) ]

    missing_water_ec_after = soil.df['water_ec'].isna()
    
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + (
            "--> Calculated by fitting (R2 = "+str(R2)+") Rhoades function in predict.water_ec.fitting_rhoades"
            if missing_water_ec_before[x] and not missing_water_ec_after[x]
            else "--> Provide water_ec, otherwise water and bulk_ec_dc_tc"
            if missing_water_ec_before[x] and missing_water_ec_after[x]
            else "")
        if missing_water_ec_before[x]
        else soil.info.water_ec[x]
        for x in range(soil.n_states)]
    

def fitting_hilhorst(soil):
    """
    Calculate missing values of soil.df.water_ec using the Hilhorst function in a fitting approach
    
    This function selects calibration data based on available bulk electrical conductivity, bulk permeability, and water permeability, and applies the Hilhorst function
      to estimate soil water electrical conductivity and an offset parameter for permeability. 
      It then performs a fitting process to optimize parameters using the objective function that minimizes the residuals between the calculated and observed bulk permeability.

    Parameters
    ----------
    soil : object
        A custom soil object containing:
        
        df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: bulk_ec, bulk_perm, water_perm, offset_perm, and water_ec.
        n_states : int
            The count of soil states.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        roundn : int
            Number of decimal places to round results.

    Returns
    -------
    None
        The function updates the `soil` object's `df` and `info` attributes with estimated values and additional information regarding the calculation.

    External Functions
    ------------------
    WaterPerm : Calculate or set missing values of soil.df.water_perm and return
    Hilhorst : Calculate the soil bulk real relative dielectric permittivity using the Hilhorst model and return
    R2_score : Calculate the coefficient of determination (R^2) of a prediction and return.

    Notes
    -----
    - The function targets soil states with known bulk electrical conductivity and bulk permeability greater than or equal to 10.
    - A least squares optimization is used to find the best parameters that fit the Hilhorst function to the calibration data.
    """
    WaterPerm(soil)

    # Selecting calibration data
    arg_EC_wn = np.array([soil.df.bulk_ec_dc_tc[x] if not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                            else np.nan for x in range(soil.n_states)])
    arg_bulk_perm_wn = np.array([soil.df.bulk_perm[x] if not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                              else np.nan for x in range(soil.n_states)])
    arg_water_perm_wn = np.array([soil.df.water_perm[x] if not np.isnan(soil.df.bulk_ec_dc_tc[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
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
    def objective_water_ec(param, bulk_perm, bulk_ec_dc_tc, water_perm):
        water_ec, offset_perm = param
        residuals = (Hilhorst(bulk_ec_dc_tc, water_ec, water_perm, offset_perm) - bulk_perm)**2
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
    
    missing_water_ec_before = soil.df['water_ec'].isna()

    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) 
                           if np.isnan(soil.df.water_ec[x]) 
                           else soil.df.water_ec[x] for x in range(soil.n_states) ]

    missing_water_ec_after = soil.df['water_ec'].isna()
    
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + (
            "--> Calculated by fitting (R2="+str(R2)+") Hilhorst function in predict.water_ec.fitting_hilhorst"
            if missing_water_ec_before[x] and not missing_water_ec_after[x]
            else "--> Provide water_ec; otherwise, bulk_perm and bulk_ec_dc_tc"
            if missing_water_ec_before[x] and missing_water_ec_after[x]
            else "")
        if missing_water_ec_before[x]
        else soil.info.water_ec[x]
        for x in range(soil.n_states)]
