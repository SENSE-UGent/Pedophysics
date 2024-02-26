import numpy as np
from scipy.optimize import minimize

from pedophysics.utils.stats import R2_score
from pedophysics.pedophysical_models.bulk_ec import Fu, WunderlichEC

from .water_ec import WaterEC
from .porosity import Porosity
from .solid_ec import SolidEC
from .frequency_ec import FrequencyEC
from .texture import Texture


def WaterFromEC(soil):
    """ 
    Calculate missing values of soil.df.water based on soil.df.bulk_ec_dc_tc 

    This function evaluates the availability of water content and bulk electrical conductivity data (bulk_ec_dc_tc) across soil states. 
    A fitting approach is applied if there are at least three soil states with known water content and bulk electrical conductivity. 
    A non-fitting approach is considered when water content is unknown and bulk electrical conductivity is known for any soil state.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A pandas DataFrame containing the soil states with columns for `water` and `bulk_ec_dc_tc`.
        - n_states: int
            The number of soil states represented in the `df`.

    Returns
    -------
    None

    Notes
    -----
    - The fitting approach requires at least three soil states with known water content and bulk electrical conductivity for reliable estimation.
    - The non-fitting approach is applied to individual soil states where water content is unknown but bulk electrical conductivity is available.

    External Functions
    ------------------
    FrequencyEC : Set missing values of soil.df.frequency_ec and return 
    fitting : Calculate missing values of soil.df.water using a fitting approach.
    non_fitting : Calculate missing values of soil.df.water using a non-fitting approach.

    Example
    -------
    >>> sample = Soil( bulk_ec = [0.01, np.nan, 0.025, 0.030, 0.040],
                clay = 10,
                porosity = 0.47,
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


def non_fitting(soil):
    """ 
    Calculate missing values of soil.df.water using a non-fitting approach.


    This function applies the Fu function within a minimization process to estimate soil water content based on soil properties such as 
    clay content, porosity, water electrical conductivity (EC), solid EC, dry EC, and saturated EC. 
    The estimation is performed for each soil state where water content is unknown.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A pandas DataFrame containing the soil states with columns for `clay`, `porosity`, `water_ec`, `solid_ec`, `dry_ec`, `sat_ec`, `bulk_ec_dc_tc`, and potentially `water`.
        - n_states: int
            The number of soil states represented in the `df`.
        - roundn: int
            The number of decimal places for rounding estimated water content values.

    bulk_ec_dc : array-like
        Soil bulk real electrical conductivity at DC frequency [S/m].

    Notes
    -----
    - The Fu function is utilized in a minimization process to estimate water content by minimizing the difference between the estimated and actual bulk ECDCTC.
    - The estimation process is applied to each soil state where water content is unknown.


    External functions
    --------
    Fu: Calculate the soil bulk real electrical conductivity using the Fu model and return
    Texture: Calculate missing values of soil.df.sand, soil.df.silt, and soil.df.clay and return
    Porosity: Calculate missing values of soil.df.porosity and return
    WaterEC: Compute missing values of soil.df.water_ec and return  
    SolidEC: Set missing values of soil.df.solid_ec and return
    """    
    Texture(soil)
    Porosity(soil)
    WaterEC(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water using Fu
    def objective_func_wat(x, clay, porosity, water_ec, solid_ec, dry_ec, sat_ec, EC):
        return (Fu(x, clay, porosity, water_ec, solid_ec, dry_ec, sat_ec) - EC)**2
    wat = []

    # Calculating water
    for i in range(soil.n_states):
        
        res = minimize(objective_func_wat, 0.15, args=(soil.df.clay[i], soil.df.porosity[i], soil.df.water_ec[i], soil.df.solid_ec[i], 
                                                        soil.df.dry_ec[i], soil.df.sat_ec[i], soil.df.bulk_ec_dc_tc[i]), bounds=[(0, .65)] )
        wat.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

   # Saving calculated water and its info
    soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_from_ec.non_fitting" if  np.isnan(soil.df.water[x]) or 
                          soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_from_ec.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
    
    soil.df['water'] = [round(wat[i], soil.roundn) if np.isnan(soil.df.water[i]) else soil.df.water[i] for i in range(soil.n_states) ]


def fitting(soil):
    """ 
    Calculate missing values of soil.df.water using a fitting approach.

    This function evaluates soil states with known water content and bulk electrical conductivity to determine initial parameters for the WunderlichEC model. 
    If the Lw parameter associated with the model is unknown, it is optimized based on the root mean square error (RMSE) between estimated and actual bulk electrical conductivity. 
    Water content is then estimated for all soil states within a valid bulk electrical conductivity range using the optimized Lw parameter and the WunderlichEC model.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A pandas DataFrame containing the soil states with columns for `water`, `bulk_ec_dc_tc`, `water_ec`, and potentially `Lw`.
        - n_states: int
            The number of soil states represented in the `df`.
        - range_ratio: float
            A ratio used to determine the range of bulk electrical conductivity values considered valid for the model.
        - roundn: int
            The number of decimal places for rounding calculated water content values.
        - Lw: float or np.nan
            The WunderlichEC model parameter, if known; otherwise, np.nan.

    Returns
    -------
    None
        The function directly modifies the `soil` object's `df` and `info` attributes with the estimated water content and does not return any value.

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    The function either estimates or uses the known Lw parameter for the WunderlichEC model and 
    fits the model to the calibration data.

    External Functions
    ------------------
    WunderlichEC: Calculate the soil bulk real electrical conductivity using the Wunderlich model and return
    WaterEC: Compute missing values of soil.df.water_ec and return  
    """

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
