import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP, LR, LR_W, LR_MV
from pedophysics.utils.stats import R2_score

from .water_perm import *
from .frequency_perm import *
from .bulk_ec import *
from .bulk_perm_inf import *
from .porosity import *
from .air_perm import *
from .solid_perm import *
from .temperature import *
from .texture import *

def BulkPerm(soil):
    """ 
    Computes missing values of soil.df.bulk_perm and return.

    Uses the given electromagnetic frequency to determine the method to estimate bulk permittivity. 
    If no frequency is provided or all frequencies are the same, different methods or conditions 
    are invoked.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_perm and frequency_perm.
        - n_states : int
            Number of soil states
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    Returns
    -------
    numpy.ndarray
        soil.df.bulk_perm.values: Array containing updated soil bulk real relative dielectric permittivity values

    Notes
    -----
    The function modifies the soil object in-place by updating the `df` and `info` attributes based 
    on the given conditions and returns the updated values of bulk permittivity.

    External Functions
    ------------------
    - FrequencyPerm : Set missing values of soil.df.frequency_perm and return 
    - fixed_freq, changing_freq : Functions invoked based on the condition of the EM frequency provided.

    Example
    -------
    >>> sample = Soil( water = [0.3, 0.1, 0.15, 0.23, 0.02],
                porosity = 0.434,
                texture = 'Silt loam',
                instrument = 'GPR')

    >>> predict.BulkPerm(sample) 
    array([18.563,  6.781,  9.328, 13.973,  3.357])
    """

    if (np.isnan(soil.df.bulk_perm)).any():  # Go over if any value is missing        
        FrequencyPerm(soil)

        # Condition to ask for frequency data
        if (np.isnan(soil.df.frequency_perm)).all():
            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Unmodified value. Please provide soil.frequency_perm" for x in range(soil.n_states)]

        # Condition for fixed EM frequency
        elif np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]):
            fixed_freq(soil)

        # Condition for changing EM frequency
        else:
            changing_freq(soil)

    return soil.df.bulk_perm.values


def fixed_freq(soil):
    """ 
    Decide between fitting and non-fitting approaches.

    Based on the given soil state data, chooses between two methodologies - fitting or non-fitting 
    to determine the bulk permittivity values for a fixed frequency.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_perm, water, and frequency_perm.
        - n_states : int
            Number of soil states

    Returns
    -------
    None
        
    Notes
    -----
    The function modifies the soil object in-place based on the given conditions.

    External Functions
    ------------------
    - fitting : Calculate missing values of soil.df.bulk_perm using a fitting approach
    - non_fitting : Calculate missing values of soil.df.bulk_perm using a non fitting approach.
    """

    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_perm[x]) for x in range(soil.n_states)) >= 3 :
        fitting(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.bulk_perm[x]) and not np.isnan(soil.df.water[x]) and soil.df.frequency_perm[x] >= 5 and soil.df.frequency_perm[x] <= 30e9 for x in range(soil.n_states) ):
        non_fitting(soil)

        
def fitting(soil):
    """ 
    Calculate missing values of soil.df.bulk_perm using a fitting approach

    This function fits the WunderlichP function to the soil data to determine the bulk permittivity values.
    It can either calculate the Lw parameter or use a provided value. The accuracy of the fitting is 
    determined by the R2 score. 

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, water_perm, and bulk_perm
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - Lw : float
            Soil scalar depolarization factor of water aggregates (effective medium theory)
        - roundn : int
            Number of decimal places to round results.
        - range_ratio : float
            Ratio to extend the domain of the regression by fitting approach.
        - n_states : int
            Number of soil states

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes and does not return any value.

    Notes
    -----
    The function modifies the soil object in-place by updating the `df` and `info` attributes based on the fitting results.
    
    External Functions
    ------------------
    - Temperature : Set missing values of soil.df.temperature and return 
    - WaterPerm : Calculate or set missing values of soil.df.water_perm and return
    - WunderlichP : Calculate the soil bulk real relative dielectric permittivity using the Wunderlich model and return
    """

    Temperature(soil)
    WaterPerm(soil)                      

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_perm) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_perm_init = min(soil.df.bulk_perm[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0

    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.bulk_perm, WunderlichP(soil.df.water, bulk_perm_init, water_init, soil.df.water_perm, soil.Lw)), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid water range
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.bulk_perm.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_perm[x]))
                                or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.bulk_perm.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
        soil.df['bulk_perm'] = [round(WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], soil.Lw), soil.roundn) 
                              if (min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]


def non_fitting(soil):
    """ 
    Calculate missing values of soil.df.bulk_perm using a non fitting approach

    This function determines the bulk permittivity of soil based on given conditions and known empirical relationships. 
    Depending on the electromagnetic (EM) frequency range provided, various functions are used.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water, bulk_perm, frequency_perm, bulk_ec, bulk_perm_inf, porosity, clay, air_perm, solid_perm, water_perm, and CEC.
        - roundn : int
            Number of decimal places to round results.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes and does not return any value.

    Notes
    -----
    The function modifies the soil object in-place by updating the `df` and `info` attributes based on the chosen approach.
    
    External Functions
    ------------------
    - BulkPermInf : Determines the bulk permittivity at infinite frequency.
    - BulkEC : Estimates the soil's bulk electrical conductivity.
    - Temperature : Set missing values of soil.df.temperature and return 
    - Porosity : Calculate missing values of soil.df.porosity and return
    - AirPerm : Set missing values of soil.df.air_perm and return
    - SolidPerm : Set missing values of soil.df.solid_perm and return
    - WaterPerm : Calculate or set missing values of soil.df.water_perm and return
    - Texture : Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    - LongmireSmithP : Calculate the soil bulk real relative dielectric permittivity using the Longmire-Smith model and return
    - LR_MV : Calculate the soil volumetric water content using the Lichtenecker and Rother model modified by Mendoza-Veirana and return
    - LR : Calculate the soil volumetric water content using the Lichtenecker and Rother model.
    - LR_W : Calculate the soil volumetric water content using the Lichtenecker and Rother model modified by Wunderlich and return
    """

    # Condition for lowest EM frequency
    if any(np.isnan(soil.df.bulk_perm[x]) and soil.df.frequency_perm[x] >= 5 and soil.df.frequency_perm[x] < 30e6 for x in range(soil.n_states)):

        BulkPermInf(soil)              
        BulkECDC(soil)

        # Saving calculated bulk_perm and its info
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.non_fitting"
                                else soil.info.bulk_perm[x] for x in range(soil.n_states)]
        
        soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_ec_dc[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) 
                                if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

    # Condition for EM frequency of common moisture sensors and GPR
    elif (np.isnan(soil.df.bulk_perm)).any() & ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm <= 30e9)).all(): 
        Temperature(soil)
        Porosity(soil)                      
        AirPerm(soil)                   
        SolidPerm(soil)                  
        WaterPerm(soil)               
        Texture(soil)                    

        if ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():
            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR_MV (reported R2=0.93) function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR_MV (reported R2=0.93) function in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR_MV(soil.df.water[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

        elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

            if np.isnan(soil.alpha): soil.alpha = 0.5 
            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR(soil.df.water[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

        elif ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 30e9)).all(): 

            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR_W function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR_W function in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR_W(soil.df.water[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)] 


def changing_freq(soil):
    """ 
    calculate missing values of soil.df.bulk_perm based on soil.df.bulk_ec_dc

    This function determines the bulk permittivity of soil based on given conditions and the LongmireSmithP pedophysical model. 
    It specifically addresses cases where the electromagnetic (EM) frequency is not constant.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_perm, frequency_perm, bulk_ec_dc, clay, sand, and bulk_perm_inf.
        - roundn : int
            Number of decimal places to round results.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes and does not return any value.
        
    Notes
    -----
    The function modifies the soil object in-place by updating the `df` and `info` attributes.

    External Functions
    ------------------
    - BulkPermInf : Set missing values of soil.df.bulk_perm_inf and return
    - BulkECDC : Compute missing values of soil.df.bulk_ec_dc and return
    - LongmireSmithP : Calculate the soil bulk real relative dielectric permittivity using the Wunderlich model and return
    """
    Texture(soil)
    BulkPermInf(soil)             
    BulkECDC(soil)

    # Warn about applying LongmireSmithP function to non-validated soil conditions
    if any( soil.df.frequency_perm[x] > 200e6 and (soil.df.clay[x] > 10 or soil.df.sand[x] < 90) and np.isnan(soil.df.bulk_perm[x]) for x in range(soil.n_states) ):
        states_warns = []
        for x in range(soil.n_states):
            if (soil.df.frequency_perm[x] > 200e6 and (soil.df.clay[x] > 10 or soil.df.sand[x] < 90) and np.isnan(soil.df.bulk_perm[x])):
                states_warns.append(x) 

        warnings.warn(f"LongmireSmithP function is applied to soil states {states_warns} with conditions soil.df.frequency_perm > 200e6 and soil.df.clay > 10 or (soil.df.sand < 90), for which the validation of such model is uncertain. ")

    # Saving calculated bulk_perm and its info
    soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.changing_freq" if np.isnan(soil.df.bulk_perm[x])
                            or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.changing_freq"
                            else soil.info.bulk_perm[x] for x in range(soil.n_states)]
    
    soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_ec_dc[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) 
                            if (np.isnan(soil.df.bulk_perm[x])) else soil.df.bulk_perm[x] for x in range(soil.n_states)]