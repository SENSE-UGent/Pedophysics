import numpy as np
from scipy.optimize import minimize

from .frequency_ec import *
from .particle_density import *
from .solid_ec import *
from .temperature import *
from .bulk_ec_dc import non_dc_to_dc

from pedophysics.pedophysical_models.bulk_ec import Fu, SheetsHendrickx, WunderlichEC
from pedophysics.utils.stats import R2_score


def BulkECDCTC(soil):
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


def non_dc_non_tc_to_dc_tc(soil):
    """
    
    """
    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec in predict.bulk_ec_dc_tc.non_dc_non_tc_to_dc_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec in predict.bulk_ec_dc_tc.non_dc_non_tc_to_dc_tc"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec_dc_tc'] = [soil.df.bulk_ec[x] if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    

def non_tc_to_tc(soil):
    """
    
    """
    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec_dc_tc.non_tc_to_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec_dc_tc.non_tc_to_tc"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
     
    soil.df['bulk_ec_dc_tc'] = [soil.df.bulk_ec_dc[x] if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] == 298.15 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]

    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using SheetsHendrickx function in predict.bulk_ec_dc_tc.non_tc_to_tc" if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] != 298.15
                        or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using SheetsHendrickx function in predict.bulk_ec_dc_tc.non_tc_to_tc" else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec_dc_tc'] = [SheetsHendrickx(soil.df.bulk_ec_dc[x], soil.df.temperature[x]) if np.isnan(soil.df.bulk_ec_dc_tc[x]) and soil.df.temperature[x] != 298.15 else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]


def shift_to_bulk_ec_dc_tc(soil):
    """
    
    """    
    if any(((not np.isnan(soil.df.bulk_ec[x])) or (not np.isnan(soil.df.bulk_ec_dc[x]))) and np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
        non_dc_to_dc(soil)
        non_dc_non_tc_to_dc_tc(soil)
        non_tc_to_tc(soil)


#def dc_freq(soil):
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
            Includes: water, frequency_ec, bulk_ec, bulk_dc_ec
        - n_states : int
            Number of soil states.

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
#    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.df.bulk_ec_dc[x]) for x in range(soil.n_states))>= 3:
#        fitting(soil)

    # Condition for non-fitting routine 
#    if any(not np.isnan(soil.df.water[x]) and np.isnan(soil.df.bulk_ec_dc[x])  for x in range(soil.n_states)):
#        non_fitting(soil)


def fitting(soil):
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
            Includes: water, water_ec, and bulk_ec_dc
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
    Return and compute bulk_ec_dc using a non-fitting approach.

    This function employs the Fu function (reported with an R^2 of 0.98) to estimate the 
    soil's bulk real electrical conductivity at DC frequency based on volumetric water content.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, clay, bulk_density, particle_density, bulk_ec, water_ec, solid_ec, dry_ec, sat_ec, and bulk_ec_dc.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

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
    from .water_ec import WaterEC # Lazy import to avoid circular dependency

    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    soil.info['bulk_ec_dc_tc'] = [str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec_dc_tc.non_fitting" if np.isnan(soil.df.bulk_ec_dc_tc[x]) 
                            or soil.info.bulk_ec_dc_tc[x] == str(soil.info.bulk_ec_dc_tc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec_dc_tc.non_fitting"
                            else soil.info.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
 
    soil.df['bulk_ec_dc_tc'] = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) 
                             if np.isnan(soil.df.bulk_ec_dc_tc[x]) else soil.df.bulk_ec_dc_tc[x] for x in range(soil.n_states)]
