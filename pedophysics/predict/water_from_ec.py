import numpy as np
from scipy.optimize import minimize

from pedophysics.utils.stats import R2_score
from pedophysics.pedophysical_models.bulk_ec import Fu, WunderlichEC

from .water_ec import WaterEC
from .particle_density import ParticleDensity
from .solid_ec import SolidEC
from .frequency_ec import FrequencyEC
from .texture import Texture


def WaterFromEC(soil):
    """ 
    Compute missing values of soil.df.water based on soil.df.bulk_ec.

    The function determines if the provided electrical conductivity is at non-DC or DC frequency. 
    If non-DC, it converts it to DC frequency. Then, it predicts the soil water content based on this 
    DC frequency electrical conductivity.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water, bulk_ec, and frequency_ec.
        - n_states : int
            Number of states or records in the dataframe.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.

    External functions
    ------
    FrequencyEC: Function to compute soil.df.frequency_ec missing values

    Example
    -------
    >>> sample = Soil( bulk_ec = [0.01, np.nan, 0.025, 0.030, 0.040],
                clay = 10,
                bulk_density = 1.4,
                water_ec = 0.5)

    >>> WaterFromEC(sample) 
    >>> sample.df.water
    0    0.105
    1    Nan
    2    0.185
    3    0.206
    4    0.243
    Name: water, dtype: float64
    """
    FrequencyEC(soil)

    # Check for conditions to use a fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)) >= 3:
        fitting(soil)

    # Check for conditions to use a non-fitting approach
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
        non_fitting(soil)


#def dc_freq(soil):
    """ 
    Decide between fitting and non-fitting approaches to calculate soil.df.water

    Based on the frequency of the electrical conductivity measuments, this function determines 
    whether to employ a fitting or non-fitting approach to estimate the soil's 
    volumetric water content.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state.
            Includes: water
        - n_states : int
            Number of states or records in the dataframe.

    bulk_ec_dc : array-like
        Array containing soil bulk real electrical conductivity at DC frequency [S/m] for each soil state.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function checks the number of non-missing records for water and bulk_ec_dc. If there are at least three,
    the fitting approach is taken, otherwise the non-fitting approach is considered.

    External functions
    --------
    fitting: Function that employs a fitting approach for prediction.
    non_fitting: Function that employs a non-fitting approach for prediction.
    """

    # Check for conditions to use a fitting approach
#    if sum(not np.isnan(soil.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states)) >= 3:
#        fitting(soil, bulk_ec_dc)

    # Check for conditions to use a non-fitting approach
#    if any(np.isnan(soil.df.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states)):
#        non_fitting(soil, bulk_ec_dc)


def non_fitting(soil):
    """ 
    Return and compute soil.df.water using a non-fitting approach.

    This function employs the Fu function (reported with an R^2 of 0.98) to estimate the 
    soil's volumetric water content based on its bulk real electrical conductivity at DC frequency.

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
    print('non fitting')
    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water using Fu
    def objective_func_wat(x, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec, EC):
        return (Fu(x, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec) - EC)**2
    wat = []

    # Calculating water
    for i in range(soil.n_states):
        res = minimize(objective_func_wat, 0.15, args=(soil.df.clay[i], soil.df.bulk_density[i], soil.df.particle_density[i], soil.df.water_ec[i], soil.df.solid_ec[i], 
                                                        soil.df.dry_ec[i], soil.df.sat_ec[i], soil.df.bulk_ec_dc_tc[i]), bounds=[(0, .65)] )
        wat.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

   # Saving calculated water and its info
    soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_from_ec.non_fitting" if  np.isnan(soil.df.water[x]) or 
                          soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_from_ec.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
    
    soil.df['water'] = [round(wat[i], soil.roundn) if np.isnan(soil.df.water[i]) else soil.df.water[i] for i in range(soil.n_states) ]


def fitting(soil):
    """ 
    Computes soil.df.water using a fitting approach.

    This function utilizes the WunderlichEC model to estimate the soil's volumetric water 
    content based on its electrical conductivity at DC frequency. It calculates the model's 
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
    print('fitting')

    WaterEC(soil) 
    
    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_ec_dc_tc) # States where calibration data are
    water_init = np.nanmin(soil.df.water[valids])
    bulk_ec_init = np.nanmin(soil.df.bulk_ec_dc_tc[valids])
    bulk_ec_final = np.nanmax(soil.df.bulk_ec_dc_tc[valids])
    bulk_ec_range = [round(bulk_ec_init - (bulk_ec_final-bulk_ec_init)/soil.range_ratio, soil.roundn), 
                     round(bulk_ec_final + (bulk_ec_final-bulk_ec_init)/soil.range_ratio, soil.roundn)]
    if bulk_ec_range[0] < 0:
        bulk_ec_range[0] = 0

    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain water
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_ec_dc_tc)**2))
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
            Wat_RMSE = np.sqrt((WunderlichEC(wat, bulk_ec_init, water_init, soil.df.water_ec[i], soil.Lw) - soil.df.bulk_ec_dc_tc[i])**2)
            return Wat_RMSE
        
        # Looping over soil states to obtain water using WunderlichEC function
        for i in range(soil.n_states):
            if (min(bulk_ec_range) <= soil.df.bulk_ec_dc_tc[i] <= max(bulk_ec_range)) & ~np.isnan(soil.df.bulk_ec_dc_tc[i]):
                result = minimize(objective_wat, 0.15, args=(i), bounds=[(0, .65)], method='L-BFGS-B')
                Wat_wund.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn))

            else:
                Wat_wund.append(np.nan)

        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.water[valids], np.array(Wat_wund)[valids]), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid bulk_ec range
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.water_from_ec.fitting, for soil.bulk_ec values between: "+str(bulk_ec_range) 
                              if min(bulk_ec_range) <= soil.df.bulk_ec_dc_tc[x] <= max(bulk_ec_range) and np.isnan(soil.df.water[x])
                                or soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.water_from_ec.fitting, for soil.bulk_ec values between: "+str(bulk_ec_range)
                                else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [Wat_wund[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Converting negative results due to fitting to zero
    soil.info['water'] = [str(soil.info.water[x]) + "--> Set to 0 because of < 0 results" if soil.df.water[x]<0 or soil.info.water[x] == str(soil.info.water[x]) + "--> Set to 0 because of < 0 results"
                            else soil.info.water[x] for x in range(soil.n_states)]
    
    soil.df['water'] = [ 0 if soil.df.water[x]<0 else soil.df.water[x] for x in range(soil.n_states)] 
