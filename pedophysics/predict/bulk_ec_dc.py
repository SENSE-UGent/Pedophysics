import numpy as np
from scipy.optimize import minimize

from .water_ec import *
from .frequency_ec import *
from .particle_density import *
from .solid_ec import *
from .bulk_ec_tc import *
from .bulk_ec_dc import *
from .bulk_ec_dc_tc import *

from pedophysics.pedophysical_models.bulk_ec import WunderlichEC, LongmireSmithEC, Fu


def DC_from_bulk_ec(soil):
    """
    
    """        
    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> Assumed equal to soil.df.bulk_ec because soil.df.frequency_ec[x] >= 5" if any(np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) and soil.df.frequency_ec[x] >= 5) or 
            soil.info.bulk_ec_dc[x] == str(soil.info.bulk_ec_dc[x]) + "--> Assumed equal to soil.df.bulk_ec because soil.df.frequency_ec[x] >= 5" else soil.info.bulk_ec_dc[x] for x in range(soil.n_states)]

    soil.df['bulk_ec_dc'] = soil.df.bulk_ec


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
            Includes: frequency_ec, bulk_ec and bulk_dc_ec
        - n_states : int
            Number of soil states.
        - roundn : int
            Number of decimal places to round results.
        - info : dict
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

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
    bulk_ec_dc = []

    for i in range(soil.n_states):
        if soil.df.frequency_ec[i] <= 5 or np.isnan(soil.df.bulk_ec[i]):
            bulk_ec_dc.append(soil.df.bulk_ec[i])

        elif soil.df.frequency_ec[i] > 5 and not np.isnan(soil.df.bulk_ec[i]):
            res = minimize(objective_func_ec_dc, 0.05, args=(soil.df.frequency_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 1)])
            bulk_ec_dc.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2) )

    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" if 
                    soil.df.frequency_ec[x] > 5 and not np.isnan(soil.df.bulk_ec[x]) and np.isnan(soil.df.bulk_ec_dc[x]) or 
                    soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" else soil.info.bulk_ec[x] for x in range(soil.n_states)]
            
    soil.df['bulk_ec_dc'] = [bulk_ec_dc[x] if np.isnan(soil.df.bulk_ec_dc[x]) else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)] 


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

    WaterEC(soil)                    

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_ec_dc) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_ec_init = min(soil.df.bulk_ec_dc[valids])
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
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_ec_dc)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.bulk_ec_dc, WunderlichEC(soil.df.water, bulk_ec_init, water_init, soil.df.water_ec, soil.Lw)), soil.roundn)

        # Saving calculated bulk_ec and its info with R2 and valid water range
        soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_ec_dc[x]))
                                or soil.info.bulk_ec_dc[x] == str(soil.info.bulk_ec_dc[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_ec_dc[x] for x in range(soil.n_states)]
                
        soil.df['bulk_ec_dc'] = [round(WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], soil.Lw), soil.roundn+3) if 
                      (min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_ec_dc[x]) else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)]
       

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

    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting" if np.isnan(soil.df.bulk_ec_dc[x]) 
                            or soil.info.bulk_ec_dc[x] == str(soil.info.bulk_ec_dc[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting"
                            else soil.info.bulk_ec_dc[x] for x in range(soil.n_states)]
 
    soil.df['bulk_ec_dc'] = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) 
                             if np.isnan(soil.df.bulk_ec_dc[x]) else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)]

