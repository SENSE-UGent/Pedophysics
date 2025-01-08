import numpy as np
from scipy.optimize import minimize
import warnings

from pedophysics.utils.stats import R2_score
from pedophysics.pedophysical_models.water import LR, LR_W, LR_MV
from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP

from .bulk_perm_inf import BulkPermInf
from .porosity import Porosity
from .air_perm import AirPerm
from .solid_perm import SolidPerm
from .water_perm import WaterPerm
from .texture import Texture


def WaterFromPerm(soil):
    """ 
    Calculate missing values of soil.df.water based on soil.df.bulk_perm

    This function checks if the permittivity frequency (`frequency_perm`) across all soil states is constant. 
    If it is, a specific adjustment is made using the `fixed_freq` function. 
    If the frequencies vary across soil states, a different adjustment is applied using the `changing_freq` function.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water and frequency_perm.
        - info : DataFrame
            Data Frame containing descriptive information about how each attribute was determined or modified.
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    None

    External functions
    --------
    fixed_freq: Decide between fitting and non-fitting approaches to calculate soil.df.water
    changing_freq: Calculate missing values of soil.df.bulk_dc_ec when soil.df.frequency_perm is not constant

    Example
    -------
    >>> sample = Soil(frequency_perm = 1e9, 
                    clay = 15,             
                    porosity = 0.434,
                    bulk_perm = [8, 10, 15])
    >>> WaterFromPerm(sample) 
    >>> sample.df.water
    0    0.125
    1    0.162
    2    0.246
    Name: water, dtype: float64
    """
    # Condition for constant permittivity frequency
    if np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]):
        fixed_freq(soil)

    # Condition for changing permittivity frequency
    else:
        changing_freq(soil)


def changing_freq(soil):    
    """ 
    Calculate missing values of soil.df.bulk_dc_ec when soil.df.frequency_perm is not constant.

    This function iterates through soil states to calculate the bulk EC for states where it is undefined. 
    The calculation is performed by minimizing the difference between the calculated permeability using the Longmire-Smith P function and the known bulk permeability. 
    Warnings are issued for soil states where the Longmire-Smith P function's applicability is uncertain due to soil conditions.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: frequency_perm, frequency_ec, bulk_ec_dc, and bulk_perm.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of soil states.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes and does not return any value.

    External Functions
    ------------------
    Texture : Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    BulkPermInf : Set missing values of soil.df.bulk_perm_inf and return
    LongmireSmithP : Calculate the soil bulk real relative dielectric permittivity using the Wunderlich model and return

    Notes
    -----
    - The function uses the `scipy.optimize.minimize` method for the minimization process.
    - Warnings are issued for soil states where the frequency of permittivity exceeds 200 MHz and either clay content is above 10% or sand content is below 90%, as the validity of the Longmire-Smith P model is uncertain in these conditions.


    """
    Texture(soil)
    BulkPermInf(soil)    
    bulk_ec_dc = []

    # Defining minimization function to obtain EC
    def objective(bulk_ec_dc, perm_inf, freq, bulk_perm):
        LS_perm = LongmireSmithP(bulk_ec_dc, perm_inf, freq)
        return (LS_perm - bulk_perm)**2

    # Calculating bulk EC from bulk perm when unknown
    for x in range(soil.n_states):
        if np.isnan(soil.df.bulk_ec_dc[x]):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x], soil.bulk_perm[x]), bounds=[(1e-6, 1)])
            bulk_ec_dc.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))
        else:
            bulk_ec_dc.append(np.nan)

    def warn_states(soil):
        # Warn about applying LongmireSmithP function to non-validated soil conditions
        mask_invalid = (
            (soil.df.frequency_perm > 200e6) &
            (soil.df.water > 0.22) &
            (soil.df.porosity > 0.255) &
            (soil.df.water_ec > 3.3) | (soil.df.water_ec < 0.0016) &
            (soil.df.clay > 10) | (soil.df.sand < 85) &
            np.isnan(soil.df.bulk_ec_dc)
        )

        # Find the indices of the invalid states
        states_warns = np.where(mask_invalid)[0]

        if states_warns.size > 0:
            warnings.warn(
                f"LongmireSmithP function is applied to soil states {states_warns} with conditions "
                f"frequency_perm > 200e6, and water > 0.22, and porosity > 0.255, and water_ec > 3.3 or water_ec < 0.0016, and clay > 10 or sand < 85"
                f", for which the validity of such model is uncertain."
                )   
            
    warn_states(soil)    
    missing_bulk_ec_dc_before = soil.df['bulk_ec_dc'].isna()
    soil.df['bulk_ec_dc'] = [bulk_ec_dc[x] 
                             if np.isnan(soil.df.bulk_ec_dc[x]) 
                             else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)]

    missing_bulk_ec_dc_after = soil.df['bulk_ec_dc'].isna()
    
    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + (
                                "--> Calculated using LongmireSmithP function in predict.water_from_perm.changing_freq" 
                                if missing_bulk_ec_dc_before[x] and not missing_bulk_ec_dc_after[x]
                                else "--> Provide bulk_ec_dc; otherwise, bulk_perm" 
                                if missing_bulk_ec_dc_before[x] and missing_bulk_ec_dc_after[x]
                                else "")
                            if missing_bulk_ec_dc_before[x]
                            else soil.info.bulk_ec_dc[x]
                            for x in range(soil.n_states)]


def fixed_freq(soil):
    """ 
    Decide between fitting and non-fitting approaches to calculate soil.df.water

    This function determines the calculation approach for soil water content based on the availability of data for water content and bulk permeability, 
    as well as the range of frequency of permittivity. It applies a fitting approach if there are at least three soil states with known water content and bulk permeability. 
    Otherwise, it considers a non-fitting approach when water content is unknown, bulk permeability is known, and the frequency of permittivity falls within a specified range.

    Parameters
    ----------
    soil : object
    A custom soil object that contains:
        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: frequency_perm, water, and bulk_perm.
        - n_states : int
            Number of soil states.

    Returns
    -------
    None
        The function directly modifies the `soil` object based on the selected approach and does not return any value.

    Notes
    -----
    This function modifies the soil object in-place, using either the `fitting` or the `non_fitting` function
    depending on the criteria described above.

    External functions
    --------
    fitting: Calculate missing values of soil.df.water using a fitting approach.
    non_fitting: Calculate missing values of soil.df.water using a non-fitting approach.
    """

    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_perm[x]) for x in range(soil.n_states)) >= 3:
        fitting(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) and 5 <= soil.df.frequency_perm[x] <=30e9 for x in range(soil.n_states)):
        non_fitting(soil)


def fitting(soil):
    """ 
    Calculate missing values of soil.df.water using a fitting approach.

    This function utilizes the WunderlichP model to estimate the soil's volumetric water 
    content based on its bulk real relative dielectric permittivity at constant frequency. 
    It calculates the model's parameters and fits them to the provided calibration data.
    The accuracy of the fitting is determined by the R2 score. 

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, bulk_perm, and water_perm.
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
        The function directly modifies the `soil` object based on the selected approach and does not return any value.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function either estimates or uses the known Lw parameter for the WunderlichP model and 
    fits the model to the calibration data.

    External functions
    --------
    WunderlichP : Calculate the soil bulk real relative dielectric permittivity using the Wunderlich model and return
    WaterPerm : Calculate or set missing values of soil.df.water_perm and return
    R2_score : Calculate the coefficient of determination (R^2) of a prediction and return.
    """
    WaterPerm(soil)                   

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_perm) # States where calibration data are
    water_init = np.nanmin(soil.df.water[valids])
    bulk_perm_init = np.nanmin(soil.df.bulk_perm[valids])
    bulk_perm_final = np.nanmax(soil.df.bulk_perm[valids])
    bulk_perm_range = [round(bulk_perm_init - (bulk_perm_final-bulk_perm_init)/soil.range_ratio, soil.roundn), 
                       round(bulk_perm_final + (bulk_perm_final-bulk_perm_init)/soil.range_ratio, soil.roundn)]
    if bulk_perm_range[0] < 0:
        bulk_perm_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm.values)**2))
            return Lw_RMSE
        
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        Wat_wund = []

        # Defining minimization function to obtain water
        def objective_wat(wat, i):
            return (WunderlichP(wat, bulk_perm_init, water_init, soil.df.water_perm[i], soil.Lw) - soil.df.bulk_perm[i])**2
        
        # Looping over soil states to obtain water using WunderlichP function
        for i in range(soil.n_states):

            if min(bulk_perm_range) <= soil.df.bulk_perm[i] <= max(bulk_perm_range) and ~np.isnan(soil.df.bulk_perm[i]):
                result = minimize(objective_wat, 0.15, args=(i), bounds=[(0, .65)], method='L-BFGS-B')
                Wat_wund.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn))

            else:
                Wat_wund.append(np.nan)

        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.water, np.array(Wat_wund)), soil.roundn)

        missing_water_before = soil.df['water'].isna()  

        soil.df['water'] = [Wat_wund[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
        missing_water_after = soil.df['water'].isna()  
        
        soil.info['water'] = [str(soil.info.water[x]) + (
                "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.water_from_perm.fitting, for soil.bulk_perm values between: "+str(bulk_perm_range)
                if missing_water_before[x] and not missing_water_after[x]
                else "--> Provide water; otherwise, bulk_perm. Regression valid for bulk_perm values between"+str(bulk_perm_range)
                if missing_water_before[x] and missing_water_after[x]
                else "")
            if missing_water_before[x]
            else soil.info.water[x]
            for x in range(soil.n_states)]
        

def non_fitting(soil):
    """ 
    Return and compute soil.df.water using a non-fitting approach.

    This function estimates soil bulk electrical conductivity (EC) and water content by applying different models based on the EM frequency range. 
    For frequencies between 5 Hz and 30 MHz, the Longmire-Smith P function is used to calculate bulk EC. 
    For frequencies between 30 MHz and 100 MHz, 100 MHz and 200 MHz, and 200 MHz and 30 GHz, different linear regression models (LR_MV, LR, LR_W) are applied to estimate water content.


    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: water, bulk_perm, frequency_perm, and bulk_ec_dc for each soil state.
        - n_states : int
            Number of soil states.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.

    Returns
    -------
    None
        The function directly modifies the `soil` object based on the selected approach and does not return any value.

    Notes
    -----
    - The function chooses the estimation model based on the EM frequency range of the soil states.
    - For frequencies between 5 Hz and 30 MHz, bulk EC is estimated. 
    For higher frequencies, water content is estimated using different linear regression models tailored to specific frequency ranges.


    External Functions
    ------------------
    Porosity : Calculate missing values of soil.df.porosity and return
    AirPerm : Set missing values of soil.df.air_perm and return
    SolidPerm : Set missing values of soil.df.solid_perm and return
    WaterPerm : Calculate or set missing values of soil.df.water_perm and return
    Texture : Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    BulkPermInf : Set missing values of soil.df.bulk_perm_inf and return
    LongmireSmithP : Calculate the soil bulk real relative dielectric permittivity using the Longmire-Smith model and return
    LR_W : Calculate the soil volumetric water content using the Lichtenecker and Rother model modified by Wunderlich and return
    LR : Calculate the soil volumetric water content using the Lichtenecker and Rother model.
    LR_MV : Calculate the soil volumetric water content using the Lichtenecker and Rother model modified by Mendoza-Veirana and return
    """
    Porosity(soil)                     
    AirPerm(soil)                      
    SolidPerm(soil)                   
    WaterPerm(soil)              
    Texture(soil)                     

    # Condition for EM frequencies between 5 and 30e6
    if ((soil.df.frequency_perm >= 5) & (soil.df.frequency_perm < 30e6)).all():
        BulkPermInf(soil)

        bulk_ec_dc = []
        # Defining minimization function to obtain bulk_ec_dc using LongmireSmithP
        def objective(bulk_ec_dc, perm_inf, freq_perm, bulk_perm):
            LS_perm = LongmireSmithP(bulk_ec_dc, perm_inf, freq_perm)
            return (LS_perm - bulk_perm)**2
        
        # Calculating bulk_ec_ec
        for i in range(soil.n_states):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[i], soil.df.frequency_perm[i], soil.df.bulk_perm[i]), bounds=[(1e-6, 1)], method='L-BFGS-B')
            bulk_ec_dc.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))

        # Check for missing values
        missing_bulk_ec_dc_before = soil.df['bulk_ec_dc'].isna()
        soil.df['bulk_ec_dc'] = [bulk_ec_dc[x] if np.isnan(soil.df.bulk_ec_dc[x]) else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)]

        missing_bulk_ec_dc_after = soil.df['bulk_ec_dc'].isna()
        
        # Update info for calculated bulk_ec_dc
        soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + (
                "--> Calculated using LongmireSmithP function in predict.water_from_perm.non_fitting"
                if missing_bulk_ec_dc_before[x] and not missing_bulk_ec_dc_after[x]
                else "--> Provide bulk_ec_dc; otherwise, bulk_perm"
                if missing_bulk_ec_dc_before[x] and missing_bulk_ec_dc_after[x]
                else "")
            if missing_bulk_ec_dc_before[x]
            else soil.info.bulk_ec_dc[x]
            for x in range(soil.n_states)]
    

    # Condition for EM frequencies between 30e6 and 100e6
    elif ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():

        missing_water_before = soil.df['water'].isna()

        soil.df['water'] = [round(LR_MV(soil.df.bulk_perm[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) 
                            if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

        missing_water_after = soil.df['water'].isna()

        # Update info for calculated water
        soil.info['water'] = [str(soil.info.water[x]) + (
                "--> Calculated using LR_MV function (reported R2=0.93) in predict.water_from_perm.non_fitting"
                if missing_water_before[x] and not missing_water_after[x]
                else "--> Provide water; otherwise bulk_perm, porosity, and CEC"
                if missing_water_before[x] and missing_water_after[x]
                else "")
            if missing_water_before[x]
            else soil.info.water[x]
            for x in range(soil.n_states)]
    

    # Condition for EM frequencies between 100e6 and 200e6
    elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

        missing_water_before = soil.df['water'].isna()
    
        soil.df['water'] = [round(LR(soil.df.bulk_perm[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha), soil.roundn) 
                            if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

        missing_water_after = soil.df['water'].isna()

        # Update info for calculated water
        soil.info['water'] = [str(soil.info.water[x]) + (
                "--> Calculated using LR function (reported RMSE=0.032) in predict.water_from_perm.non_fitting"
                if missing_water_before[x] and not missing_water_after[x]
                else "--> Provide water; otherwise bulk_perm, and porosity"
                if missing_water_before[x] and missing_water_after[x]
                else "")
            if missing_water_before[x]
            else soil.info.water[x]
            for x in range(soil.n_states)]
        

    # Condition for EM frequencies between 200e6 and 30e9
    elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 30e9))).all():
      
        missing_water_before = soil.df['water'].isna()

        soil.df['water'] = [round(LR_W(soil.df.bulk_perm[x], soil.df.porosity[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) 
                            if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)] 

        missing_water_after = soil.df['water'].isna()

        # Update info for calculated water
        soil.info['water'] = [str(soil.info.water[x]) + (
                "--> Calculated using LR_W function in predict.water_from_perm.non_fitting"
                if missing_water_before[x] and not missing_water_after[x]
                else "--> Provide water; otherwise bulk_perm, porosity, and Clay"
                if missing_water_before[x] and missing_water_after[x]
                else "")
            if missing_water_before[x]
            else soil.info.water[x]
            for x in range(soil.n_states)]