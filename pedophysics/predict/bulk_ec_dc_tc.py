import numpy as np
from scipy.optimize import minimize

from .frequency_ec import *
from .solid_ec import *
from .temperature import *
from .bulk_ec_dc import non_dc_to_dc
from .texture import Texture
from .porosity import Porosity

from pedophysics.pedophysical_models.bulk_ec import Fu, SheetsHendrickx, WunderlichEC
from pedophysics.utils.stats import R2_score


def BulkECDCTC(soil):
    """ 
    Compute missing values of soil.df.bulk_ec_dc_tc and return

    This function initiates a series of steps to estimate missing `bulk_ec_dc_tc` values in the soil DataFrame.
    It starts by updating the soil's frequency and temperature data using `FrequencyEC` and `Temperature` functions,
    respectively, followed by a shift to `bulk_ec_dc_tc` values through `shift_to_bulk_ec_dc_tc`. If the conditions
    for a fitting routine are met, indicated by at least three states having non-NaN values for both water content
    and `bulk_ec_dc_tc`, the `fitting` function is called to further refine these estimates. Conversely, if there
    are states with non-NaN water content but missing `bulk_ec_dc_tc` values, the `non_fitting` routine is invoked
    to estimate these missing values.

    Parameters
    ----------
    soil : object
        A custom soil object containing:
        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_ec, bulk_ec_dc_tc.
        - n_states : int
            Number of soil states.

    Returns
    -------
    numpy.ndarray
        soil.df.bulk_ec_dc_tc.values: array containing updated soil bulk real electrical conductivity direct current (0 Hz) temperature corrected (298.15 K) values

    Notes
    -----
    - The function operates in-place, modifying the soil object's DataFrame directly.
    - It relies on both fitting and non-fitting routines to estimate missing `bulk_ec_dc_tc` values, depending
      on the availability and condition of the soil data.
    - The estimation process is contingent upon the presence of adequate non-NaN data points for both water content
      and `bulk_ec_dc_tc` across the soil states.

    External Functions
    ------------------
    - FrequencyEC : Set missing values of soil.df.frequency_ec and return 
    - Temperature : Set missing values of soil.df.temperature and return 
    - shift_to_bulk_ec_dc_tc : Compute missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec or soil.df.bulk_ec_dc
    - fitting : Calculate missing values of soil.df.bulk_ec_dc_tc using a fitting approach
    - non_fitting : Calculate missing values of soil.df.bulk_ec_dc_tc using a non-fitting approach
    """
    if (np.isnan(soil.df.bulk_ec_dc_tc)).any():  # Go over if any value is missing      

        FrequencyEC(soil)
        Temperature(soil)
        shift_to_bulk_ec_dc_tc(soil)

        # Condition for fitting routine 
        if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states))>= 3:
            fitting(soil)

        # Condition for non-fitting routine 
        if any(not np.isnan(soil.df.water[x]) and np.isnan(soil.df.bulk_ec_dc_tc[x])  for x in range(soil.n_states)):
            non_fitting(soil)

    return soil.df.bulk_ec_dc_tc.values


def shift_to_bulk_ec_dc_tc(soil):
    """
    Compute missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec or soil.df.bulk_ec_dc
    
    This function checks each soil state for the presence of non-NaN `bulk_ec` or `bulk_ec_dc` values, along with
    NaN `bulk_ec_dc_tc` values. When such a condition is met, it sequentially calls three functions to perform
    conversions and updates: `non_dc_to_dc` to convert non-direct current (non-DC) electrical conductivity (EC)
    values to direct current (DC) values, `non_dc_non_tc_to_dc_tc` to convert non-DC, non-temperature corrected (non-TC)
    EC values to DC, temperature-corrected (TC) values, and finally, `non_tc_to_tc` to convert non-temperature corrected
    EC values to temperature-corrected values.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A DataFrame containing quantitative information of soil attributes for each state, including
            `bulk_ec`, `bulk_ec_dc`, and `bulk_ec_dc_tc`.
        - n_states: int
            The number of soil states.

    Returns
    -------
    None
        This function does not return a value. It operates in-place, modifying the soil object's DataFrame directly.

    External Functions
    ------------------
    non_dc_to_dc : Calculate missing values of soil.df.bulk_ec_dc based on soil.df.bulk_ec
    non_dc_non_tc_to_dc_tc : Calculate missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec
    non_tc_to_tc : Calculate missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec_dc

    Notes
    -----
    - The function operates in-place, modifying the `soil.df` DataFrame directly.
    - It is designed to handle soil states that have either `bulk_ec` or `bulk_ec_dc` values available but are missing `bulk_ec_dc_tc` values.
    - This sequential update process ensures that all relevant EC values are converted and updated to their DC, TC equivalents where necessary.

    """    
    if any(((not np.isnan(soil.df.bulk_ec[x])) or (not np.isnan(soil.df.bulk_ec_dc[x]))) and np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
        non_dc_to_dc(soil) 
        non_dc_non_tc_to_dc_tc(soil) 
        non_tc_to_tc(soil)


def non_dc_non_tc_to_dc_tc(soil):
    """
    Calculate missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec

    If a `bulk_ec_dc_tc` value is missing (`NaN`) and the corresponding soil temperature is 298.15K with an electrical conductivity (EC) measurement frequency of 5Hz or less, 
    the `bulk_ec_dc_tc` value is set equal to the non-DC, non-temperature corrected `bulk_ec` value for that state. 
    Additionally, an annotation is added to `soil.info['bulk_ec_dc_tc']` indicating that the updated value is equal to the `bulk_ec` value under the specified conditions, 
    as part of the `predict.bulk_ec_dc_tc.non_dc_non_tc_to_dc_tc` process.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A DataFrame containing quantitative information of soil attributes for each state, including
            `bulk_ec`, `bulk_ec_dc_tc`, `temperature`, and `frequency_ec`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - n_states: int
            The number of soil states.

    Returns
    -------
    None
        This function does not return a value. It operates in-place, modifying the soil object's DataFrame and info dictionary directly.

    Notes
    -----
    - The function is specifically designed to handle cases where direct measurements of temperature-corrected DC bulk EC are not available but can be inferred from non-DC EC measurements under standard conditions (298.15K and ≤5Hz).
    - Annotations in `soil.info['bulk_ec_dc_tc']` aid in tracking the source and method used for estimating the updated `bulk_ec_dc_tc` values.
    """
    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec in predict.bulk_ec_dc_tc.non_dc_non_tc_to_dc_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec in predict.bulk_ec_dc_tc.non_dc_non_tc_to_dc_tc"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec_dc_tc'] = [soil.df.bulk_ec[x] if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    

def non_tc_to_tc(soil):
    """
    Calculate missing values of soil.df.bulk_ec_dc_tc based on soil.df.bulk_ec_dc

    This function iterates over the soil states to update `bulk_ec_dc_tc` values. 
    For states where `bulk_ec_dc_tc` is missing and the temperature is exactly 298.15K, the value is set equal to the non-temperature corrected `bulk_ec_dc`. 
    An annotation is added to `soil.info['bulk_ec_dc_tc']` indicating that the updated value is equivalent to `soil.df.bulk_ec_dc` under standard temperature conditions, 
    as part of the `predict.bulk_ec_dc_tc.non_tc_to_tc` process.
    If the temperature deviates from 298.15K and `bulk_ec_dc_tc` is missing, the `SheetsHendrickx` function is used to calculate the temperature-corrected value from `bulk_ec_dc`. 
    This is also annotated in `soil.info['bulk_ec_dc_tc']` to indicate that the value was calculated using the `SheetsHendrickx` function under non-standard temperature conditions.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A DataFrame containing quantitative information of soil attributes for each state, including
            `bulk_ec_dc`, `bulk_ec_dc_tc`, and `temperature`.
        - info: DataFrame
            A dictionary containing descriptive information about how each array-like attribute was calculated.
        - n_states: int
            The number of soil states.

    Returns
    -------
    None
        This function does not return a value. It operates in-place, modifying the soil object's DataFrame and info dictionary directly.

    External Functions
    ------------------
    SheetsHendrickx : Calculate the soil bulk real electrical conductivity using the Sheets-Hendricks model and return

    Notes
    -----
    - The function distinguishes between standard (298.15K) and non-standard temperatures for updating `bulk_ec_dc_tc` values.
    - Annotations in `soil.info['bulk_ec_dc_tc']` provide insight into the source of the updated values, enhancing data traceability.
    - The `SheetsHendrickx` function is utilized for temperature corrections under non-standard conditions, applying a model to estimate the temperature-corrected EC value.    
    """
    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec_dc_tc.non_tc_to_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec_dc_tc.non_tc_to_tc"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
     
    soil.df['bulk_ec_dc_tc'] = [soil.df.bulk_ec_dc[x] if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]

    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using SheetsHendrickx function in predict.bulk_ec_dc_tc.non_tc_to_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] != 298.15
                        or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using SheetsHendrickx function in predict.bulk_ec_dc_tc.non_tc_to_tc" else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec_dc_tc'] = [SheetsHendrickx(soil.df.bulk_ec_dc[x], soil.df.temperature[x]) if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] != 298.15 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]


def fitting(soil):
    """ 
    Calculate missing values of soil.df.bulk_ec_dc_tc using a fitting approach

    This function utilizes the WunderlichEC model to estimate the soil's bulk real electrical conductivity at DC frequency temperature-corrected based on water content. 
    It calculates the model's parameters and fits them to the provided calibration data. 
    The accuracy of the fitting is determined by the R2 score. 

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, water_ec, and bulk_ec_dc_tc
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - Lw : float
            Soil scalar depolarization factor of water aggregates (effective medium theory)
        - roundn : int
            Number of decimal places to round results.
        - range_ratio : float
            Ratio to extend the domain of the regression by fitting approach.
        - n_states : int
            Number of soil states. 

    Returns
    -------
    None 
        This function does not return a value. It operates in-place, modifying the soil object's DataFrame and info dictionary directly.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function either estimates or uses the known Lw parameter for the WunderlichEC model and 
    fits the model to the calibration data.

    External Functions
    ------------------
    WunderlichEC: Calculate the soil bulk real electrical conductivity using the Wunderlich model and return
    WaterEC: Calculate missing values of soil.df.water based on soil.df.bulk_ec_dc_tc 
    """
    from .water_ec import WaterEC # Lazy import to avoid circular dependency

    WaterEC(soil)                    

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_ec_dc_tc) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_ec_dc_tc_init = min(soil.df.bulk_ec_dc_tc[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_dc_tc_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_ec_dc_tc)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.bulk_ec_dc_tc, WunderlichEC(soil.df.water, bulk_ec_dc_tc_init, water_init, soil.df.water_ec, soil.Lw)), soil.roundn)

        # Saving calculated bulk_ec_dc_tc and its info with R2 and valid water range
        soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec_dc_tc.fitting, for soil.water values between"+str(water_range) if np.isnan(soil.df.bulk_ec_dc_tc[x]) and (min(water_range) <= soil.water[x] <= max(water_range))
                                or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec_dc_tc.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
                
        soil.df['bulk_ec_dc_tc'] = [round(WunderlichEC(soil.df.water[x], bulk_ec_dc_tc_init, water_init, soil.df.water_ec[x], soil.Lw), soil.roundn+3) if 
                      np.isnan(soil.df.bulk_ec_dc_tc[x]) and (min(water_range) <= soil.water[x] <= max(water_range)) else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]


def non_fitting(soil):
    """ 
    Calculate missing values of soil.df.bulk_ec_dc_tc using a non-fitting approach

    This function employs the Fu function (reported with an R^2 of 0.98) to estimate the 
    soil's bulk real electrical conductivity at DC frequency temperature-corrected based on volumetric water content.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, clay, porosity, bulk_ec, water_ec, solid_ec, dry_ec, sat_ec, and bulk_ec_dc_tc.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

    Returns
    -------
    None
        This function does not return a value. It operates in-place, modifying the soil object's DataFrame and info dictionary directly.

    External functions
    --------
    Fu: Calculate the soil bulk real electrical conductivity using the Fu model and return
    Texture: Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    Porosity: Calculate missing values of soil.df.porosity and return
    WaterEC: Compute missing values of soil.df.water_ec and return  
    SolidEC: Set missing values of soil.df.solid_ec and return
    
    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function uses optimization techniques to minimize the difference between the Fu function output 
    and the provided bulk real DC electrical conductivity temperature-corrected to determine the volumetric water content.
    """
    from .water_ec import WaterEC # Lazy import to avoid circular dependency

    Texture(soil)
    Porosity(soil)
    WaterEC(soil)
    SolidEC(soil)

    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec_dc_tc.non_fitting" if np.isnan(soil.df.bulk_ec_dc_tc[x]) 
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec_dc_tc.non_fitting"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
 
    soil.df['bulk_ec_dc_tc'] = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.porosity[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) 
                             if np.isnan(soil.df.bulk_ec_dc_tc[x]) else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
