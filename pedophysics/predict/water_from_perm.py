import numpy as np
from scipy.optimize import minimize

from pedophysics.utils.stats import R2_score
from pedophysics.pedophysical_models.water import LR, LR_W, LR_MV
from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP

from .bulk_perm_inf import BulkPermInf
from .particle_density import ParticleDensity
from .air_perm import AirPerm
from .solid_perm import SolidPerm
from .water_perm import WaterPerm
from .texture import Texture


def WaterFromPerm(soil):
    """ 
    Compute missing values of soil.df.water based on soil.df.bulk_perm

    Depending on the consistency of the permittivity frequency provided, this function decides 
    whether to apply a fixed frequency or changing frequency method to predict soil volumetric water content.

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

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.

    Two main methods are applied based on the nature of the permittivity frequency:
    - `fixed_freq`: applied when the permittivity frequency is constant across the soil states.
    - `changing_freq`: applied when the permittivity frequency varies across the soil states.

    External functions
    --------
    fixed_freq: Function to predict water content based on a constant permittivity frequency.
    changing_freq: Function to predict water content based on varying permittivity frequencies.

    Example
    -------
    >>> sample = Soil(frequency_perm = 1e9, 
                    clay = 15,             
                    bulk_density = 1.5,
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
    Predict soil attributes when permittivity frequency is changing.

    Determines bulk electrical conductivity (bulk_ec) using the LongmireSmithP function, given the changing 
    nature of permittivity frequency in the provided soil data.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: frequency_perm, frequency_ec, bulk_ec, and bulk_perm.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - n_states : int
            Number of soil states.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The minimization function `objective` computes the difference between the LongmireSmithP predicted 
    permittivity and actual permittivity to obtain the best bulk_ec for each soil state.

    External functions
    --------
    LongmireSmithP: Function used to predict soil bulk real relative dielectric permittivity given bulk_ec, perm_inf, and frequency.
    """

    BulkPermInf(soil)    
    bulk_ec = []

    # Defining minimization function to obtain EC
    def objective(bulk_ec, perm_inf, freq, bulk_perm):
        LS_perm = LongmireSmithP(bulk_ec, perm_inf, freq)
        return (LS_perm - bulk_perm)**2

    # Calculating bulk EC from bulk perm when unknown
    for x in range(soil.n_states):
        if np.isnan(soil.df.bulk_ec[x]):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x], soil.bulk_perm[x]), bounds=[(1e-6, 1)])
            bulk_ec.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))
        else:
            bulk_ec.append(np.nan)

    # Saving calculated bulk_ec and its info
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.changing_freq" if np.isnan(soil.df.bulk_ec[x])
                            or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.changing_freq"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
        
    soil.df['bulk_ec'] = [bulk_ec[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]

    # Saving calculated frequency_ec and its info
    soil.info['frequency_ec'] = [str(soil.info.frequency_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.changing_freq" if not np.isnan(bulk_ec[x])
                            or soil.info.frequency_ec[x] == str(soil.info.frequency_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.changing_freq"
                            else soil.info.frequency_ec[x] for x in range(soil.n_states)]
        
    soil.df['frequency_ec'] = [0 if not np.isnan(bulk_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]


def fixed_freq(soil):
    """ 
    Decide between fitting and non-fitting approaches to calculate soil.df.water

    Determines the soil's water content based on the bulk permittivity and the permittivity frequency.
    The approach to estimate water content depends on the available data: 
    1) Fitting approach: Used if there are at least 3 non-NaN values of water and bulk permittivity (calibration data).
    2) Non-fitting approach: Used if there's any soil state with NaN water, non-NaN bulk permittivity, 
       and a frequency_perm value between 5 and 30e9.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water, bulk_perm, and frequency_perm.
        - n_states : int
            Number of soil states

    Notes
    -----
    This function modifies the soil object in-place, using either the `fitting` or the `non_fitting` function
    depending on the criteria described above.

    External functions
    --------
    fitting: Function used to fit soil data and predict missing water content values.
    non_fitting: Function used to predict water content without a fitting approach.
    """

    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_perm[x]) for x in range(soil.n_states)) >= 3:
        fitting(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) and 5 <= soil.df.frequency_perm[x] <=30e9 for x in range(soil.n_states)):
        non_fitting(soil)


def fitting(soil):
    """ 
    Computes soil.df.water using a fitting approach.

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

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function either estimates or uses the known Lw parameter for the WunderlichP model and 
    fits the model to the calibration data.

    External functions
    --------
    WunderlichP: Function that defines the relationship between water content and relative dielectric permittivity.
    WaterPerm: Function to compute soil water real relative dielectric permittivity.
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
        R2 = round(R2_score(soil.df.water.values, Wat_wund), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid bulk_perm range
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.water_from_perm.fitting, for soil.bulk_perm values between: "+str(bulk_perm_range) if min(bulk_perm_range) <= soil.df.bulk_perm[x] <= max(bulk_perm_range) and ~np.isnan(soil.df.bulk_perm[x]) and np.isnan(soil.df.water[x])
                                or soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.water_from_perm.fitting, for soil.bulk_perm values between: "+str(bulk_perm_range)
                                else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [Wat_wund[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]


def non_fitting(soil):
    """ 
    Return and compute soil.df.water using a non-fitting approach.

    Uses various methods to calculate water content and bulk electrical conductivity (bulk_ec) based on 
    different electromagnetic (EM) frequency ranges and other given soil attributes.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: water, bulk_perm, frequency_perm, and bulk_ec for each soil state.
        - n_states : int
            Number of soil states.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.

    Notes
    -----
    The function works based on electromagnetic (EM) frequency conditions and modifies the soil 
    object in-place by updating the `df` and `info` attributes. Various prediction methods like LongmireSmithP,
    LR_MV, LR, and LR_W are utilized based on the frequency ranges.

    External Functions
    ------------------
    ParticleDensity, AirPerm, SolidPerm, WaterPerm, Texture, BulkPermInf : Functions used to estimate various soil attributes.
    LongmireSmithP, LR_MV, LR, LR_W : Pedophysical models for predictions based on frequency ranges and other attributes.

    """
    ParticleDensity(soil)                     
    AirPerm(soil)                      
    SolidPerm(soil)                   
    WaterPerm(soil)              
    Texture(soil)                     
    #CEC(soil)                      

    # Condition for EM frequencies between 5 and 30e6
    if ((soil.df.frequency_perm >= 5) & (soil.df.frequency_perm < 30e6)).all():
        BulkPermInf(soil)

        bulk_ec = []
        # Defining minimization function to obtain bulk_ec using LongmireSmithP
        def objective(bulk_ec, perm_inf, freq_perm, bulk_perm):
            LS_perm = LongmireSmithP(bulk_ec, perm_inf, freq_perm)
            return (LS_perm - bulk_perm)**2
        
        # Calculating bulk_ec
        for i in range(soil.n_states):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[i], soil.df.frequency_perm[i], soil.df.bulk_perm[i]), bounds=[(1e-6, 1)], method='L-BFGS-B')
            bulk_ec.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))

        # Saving calculated bulk_ec and its info
        soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.non_fitting" if np.isnan(soil.df.bulk_ec[x]) 
                                or soil.info.bulk_ec[x] ==str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water_from_perm.non_fitting" else soil.info.bulk_ec[x] for x in range(soil.n_states)]
        
        soil.df['bulk_ec'] = [bulk_ec[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 30e6 and 100e6
    elif ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR_MV function (reported R2=0.93) in predict.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                              or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR_MV function (reported R2=0.93) in predict.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]

        soil.df['water'] = [round(LR_MV(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) 
                            if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 100e6 and 200e6
    elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                        or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
                
        print('LR', [LR(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha) for x in range(soil.n_states)])        
        soil.df['water'] = [round(LR(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 200e6 and 30e9
    elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 30e9))).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR_W function in predict.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                        or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR_W function in predict.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [round(LR_W(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)] 
