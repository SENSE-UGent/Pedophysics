import numpy as np
from scipy.optimize import minimize

from .water_ec import *
from .frequency_ec import *
from .particle_density import *
from .solid_ec import *

from pedophysics.pedophysical_models.bulk_ec import WunderlichEC, LongmireSmithEC, Fu


def BulkEC(soil):
    """ 
    Return and decides for different approaches to calculate soil.df.bulk_ec.

    This function computes the bulk EC of the soil based on the provided conditions and data in the soil object.
    It covers cases where the EC is provided at different frequencies and converts these values as necessary.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: bulk_ec, frequency_ec, and water. 
        - n_states : int
            Number of soil states.

    Returns
    -------
    numpy.ndarray
        Array containing updated soil bulk real electrical conductivity values

    Notes
    -----
    The function modifies the soil object in-place by updating the `df` attribute based on the provided conditions.

    External Functions
    ------------------
    - FrequencyEC : Function to determine the EC at specific frequencies.
    - non_dc_to_dc : Converts non-direct current (DC) values to DC values for EC.
    - dc_freq : Adjusts the DC values of EC based on frequency.
    - dc_to_non_dc : Converts DC values of EC back to their non-DC counterparts.
    """

    if (np.isnan(soil.df.bulk_ec)).any():  # Go over if any value is missing        
        FrequencyEC(soil)

        if any(soil.df.frequency_ec[x] >= 5 and np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) for x in range(soil.n_states)):
            bulk_ec_dc = non_dc_to_dc(soil)

        else:
            bulk_ec_dc = soil.df.bulk_ec

        bulk_ec_dc = dc_freq(soil, bulk_ec_dc)
        dc_to_non_dc(soil, bulk_ec_dc)

    return soil.df.bulk_ec.values


def non_dc_to_dc(soil):
    """
    Convert non-direct current (non-DC) values of soil.df.bulk_ec to direct current (DC) values.

    Given the bulk EC values at various electromagnetic frequencies, this function uses the pedophysical model
    LongmireSmithEC to estimate the bulk EC of the soil at zero Hertz (direct current).

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: frequency_ec, and bulk_ec.
        - n_states : int
            Number of soil states.
        - roundn : int
            Number of decimal places to round results.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    Returns
    -------
    numpy.ndarray
        Array containing the estimated direct current (DC) values of bulk real electrical conductivity.

    Notes
    -----
    The function modifies the soil object in-place by updating the `info` attribute to include details about the 
    conversion method used for each state.

    External Functions
    ------------------
    - LongmireSmithEC : Model used to relate non-DC values to DC values.
    """

    # Defining minimization function to obtain DC bulk EC 
    def objective_func_ec_dc(bulk_ec_dc, frequency_ec, bulk_ec):
        return (LongmireSmithEC(bulk_ec_dc, frequency_ec) - bulk_ec)**2
    ec_dc = []

    for i in range(soil.n_states):
        if soil.df.frequency_ec[i] <= 5 or np.isnan(soil.df.bulk_ec[i]):
            ec_dc.append(soil.df.bulk_ec[i])

        elif soil.df.frequency_ec[i] > 5 and not np.isnan(soil.df.bulk_ec[i]):
            res = minimize(objective_func_ec_dc, 0.05, args=(soil.df.frequency_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 1)])
            ec_dc.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2) )

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" if 
                    soil.df.frequency_ec[x] > 5 and not np.isnan(soil.df.bulk_ec[x]) or 
                    soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" 
                    else soil.info.bulk_ec[x] for x in range(soil.n_states)]
            
    return np.array(ec_dc)


def dc_freq(soil, bulk_ec_dc):
    """
    Decide between fitting and non-fitting approaches to calculate soil.df.bulk_ec.

    Based on the frequency of the electrical conductivity measuments, this function determines 
    whether to employ a fitting or non-fitting approach to estimate the soil's bulk real EC.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: water, frequency_ec, and bulk_ec.
        - n_states : int
            Number of soil states.

    bulk_ec_dc : numpy.ndarray
        Array containing soil bulk real electrical conductivity at DC frequency [S/m] for each soil state.

    Returns
    -------
    numpy.ndarray
        Array containing updated soil bulk real electrical conductivity at DC frequency [S/m] for each soil state.

    Notes
    -----
    The function decides the appropriate method (fitting or non-fitting) for adjusting the bulk EC values based on
    the availability of the soil's water content data.

    External Functions
    ------------------
    - fitting : Function used to adjust bulk EC values using a fitting routine.
    - non_fitting : Function used to adjust bulk EC values using a non-fitting routine.
    """
    # Condition for fitting routine 
    if sum(not np.isnan(soil.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states))>= 3:
        bulk_ec_dc = fitting(soil, bulk_ec_dc)

    # Condition for non-fitting routine 
    if any(not np.isnan(soil.df.water[x]) and np.isnan(bulk_ec_dc[x])  for x in range(soil.n_states)):
        bulk_ec_dc = non_fitting(soil, bulk_ec_dc)

    return bulk_ec_dc


def fitting(soil, bulk_ec_dc):
    """ 
    Return and computes bulk_ec_dc using a fitting approach.

    This function utilizes the WunderlichEC model to estimate the soil's bulk real 
    electrical conductivity at DC frequency based on water content. It calculates the model's 
    parameters and fits them to the provided calibration data. The accuracy of the fitting 
    is determined by the R2 score. 

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water and water_ec.
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

    bulk_ec_dc : array-like
        Soil bulk real electrical conductivity at DC frequency [S/m].

    Returns
    -------
    numpy.ndarray
        Array containing updated soil bulk real electrical conductivity at DC frequency [S/m] for each soil state.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function either estimates or uses the known Lw parameter for the WunderlichEC model and 
    fits the model to the calibration data.

    External Functions
    ------------------
    WunderlichEC: Function that defines the relationship between water content and electrical conductivity.
    WaterEC: Function to compute soil water real electrical conductivity.
    """

    WaterEC(soil)                    

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(bulk_ec_dc) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_ec_init = min(bulk_ec_dc[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - bulk_ec_dc)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(bulk_ec_dc, WunderlichEC(soil.df.water, bulk_ec_init, water_init, soil.df.water_ec, soil.Lw)), soil.roundn)

        # Saving calculated bulk_ec and its info with R2 and valid water range
        soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_ec[x]))
                                or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
                
        bulk_ec_dc = [round(WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], soil.Lw), soil.roundn+3) if 
                      (min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(bulk_ec_dc[x]) else bulk_ec_dc[x] for x in range(soil.n_states)]
       
        return bulk_ec_dc


def non_fitting(soil, bulk_ec_dc):
    """ 
    Return and compute bulk_ec_dc using a non-fitting approach.

    This function employs the Fu function (reported with an R^2 of 0.98) to estimate the 
    soil's bulk real electrical conductivity at DC frequency based on volumetric water content.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, clay, bulk_density, particle_density, bulk_ec, water_ec, solid_ec, dry_ec, and sat_ec.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

    bulk_ec_dc : array-like
        Soil bulk real electrical conductivity at DC frequency [S/m].

    Returns
    -------
    numpy.ndarray
        Array containing updated soil bulk real electrical conductivity at DC frequency [S/m] for each soil state.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function uses optimization techniques to minimize the difference between the Fu function output 
    and the provided bulk real DC electrical conductivity to determine the volumetric water content.

    External functions
    --------
    Fu: Function that defines the relationship between water content and electrical conductivity.
    Texture: Function to derive soil texture properties.
    ParticleDensity: Function to compute particle_density.
    WaterEC: Function to compute water_ec
    SolidEC: Function to compute solid_ec
    """

    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting" if np.isnan(bulk_ec_dc[x]) 
                            or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
 
    bulk_ec_dc = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) if np.isnan(bulk_ec_dc[x]) else bulk_ec_dc[x] for x in range(soil.n_states)]

    return bulk_ec_dc


def dc_to_non_dc(soil, bulk_ec_dc):
    """
    Converts direct current (DC) bulk electrical conductivity (EC) values to non-DC frequencies.

    This function uses the LongmireSmithEC pedophysical model to adjust the direct current (DC) bulk EC values
    of the soil to the actual electromagnetic (EM) frequency. This is particularly useful when the
    actual frequency is above 5 Hz.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state.
            Includes: frequency_ec, bulk_ec, and other relevant attributes.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

    bulk_ec_dc : array-like
        Soil bulk real electrical conductivity at DC frequency [S/m].

    Notes
    -----
    The function differentiates between cases where the bulk EC value is provided by the user or calculated 
    using the LongmireSmithEC function. If the user has provided the value, it sets the 'info' attribute 
    accordingly.

    External Functions
    ------------------
    - LongmireSmithEC : Function used to adjust bulk EC values from DC to non-DC frequencies.
    """
    
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" if (np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or 
                soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                else soil.info.bulk_ec[x] for x in range(soil.n_states)]

    bulk_ec_non_dc = [round(LongmireSmithEC(bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else bulk_ec_dc[x] for x in range(soil.n_states)]

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" if (not np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or 
                soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" 
                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df["bulk_ec"] = [bulk_ec_non_dc[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]