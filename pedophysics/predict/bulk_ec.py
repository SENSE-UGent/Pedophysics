import numpy as np
from scipy.optimize import minimize

from .water_ec import *
from .frequency_ec import *
from .particle_density import *
from .solid_ec import *
from .bulk_ec_tc import *
from .bulk_ec_dc import *
from .bulk_ec_dc_tc import *

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC, SheetsHendrickx


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

        if any(np.isnan(soil.df.bulk_ec[x]) and (not np.isnan(soil.df.bulk_ec_tc[x])) for x in range(soil.n_states)):
            temperature_correction(soil)

        if any(np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) and soil.df.frequency_ec[x] >= 5 for x in range(soil.n_states)):
            non_dc_to_dc(soil)

        else:
            DC_from_bulk_ec(soil)

        dc_freq(soil)
        dc_to_non_dc(soil)

    BulkECTC(soil)
    return soil.df.bulk_ec_tc.values


def temperature_correction(soil):
    """

    """    
    # Defining minimization function to obtain DC bulk EC 
    def objective_func_ec_tc(bulk_ec, bulk_ec_tc, temperature):
        return (SheetsHendrickx(bulk_ec, temperature) - bulk_ec_tc)**2
    bulk_ec = []

    for i in range(soil.n_states):
        res = minimize(objective_func_ec_tc, 0.05, args=(soil.df.bulk_ec_tc[i], soil.df.frequency_ec[i]), bounds=[(0, 1)])
        bulk_ec.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2) )

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated from soil.df.bulk_ec_tc using SheetsHendrickx function in predict.bulk_ec.temperature_correction" if 
                    np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_ec_tc[x]) or 
                    soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated from soil.df.bulk_ec_tc using SheetsHendrickx function in predict.bulk_ec.temperature_correction"  
                    else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec'] = [round(bulk_ec[i], soil.roundn+3) if np.isnan(soil.df.bulk_ec[i]) and not np.isnan(soil.df.bulk_ec_tc[i]) else soil.df.bulk_ec[i] for i in range(soil.n_states)]


def dc_freq(soil):
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
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.df.bulk_ec_dc[x]) for x in range(soil.n_states))>= 3:
        fitting(soil)

    # Condition for non-fitting routine 
    if any(not np.isnan(soil.df.water[x]) and np.isnan(soil.df.bulk_ec_dc[x])  for x in range(soil.n_states)):
        non_fitting(soil)


def dc_to_non_dc(soil):
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
            Includes: frequency_ec, bulk_ec, bulk_ec_dc and other relevant attributes.
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

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

    bulk_ec_non_dc = [round(LongmireSmithEC(soil.df.bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else soil.df.bulk_ec_dc[x] for x in range(soil.n_states)]

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" if (not np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or 
                soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" 
                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df["bulk_ec"] = [bulk_ec_non_dc[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]