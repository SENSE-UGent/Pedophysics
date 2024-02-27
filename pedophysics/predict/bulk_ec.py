import numpy as np

from .bulk_ec_dc import BulkECDC

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC

def BulkEC(soil):
    """ 
    Calculate missing values of soil.df.bulk_ec and return

    This function checks for NaN values in the soil's bulk EC data. If any are found,
    it performs a series of operations to estimate the bulk EC, including direct current (DC)
    conversion and adjustment for non-DC components. The operations involve calling the BulkECDC function
    for DC conversion, followed by additional conversion functions to account for non-DC factors.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A pandas DataFrame containing the soil states with columns for `bulk_ec`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.

    Returns
    -------
    np.ndarray
        soil.df.bulk_ec.values: an array of updated or original soil bulk real electrical conductivity values

    External Functions
    ------------------
    BulkECDC : Compute missing values of soil.df.bulk_ec_dc and return
    conversion : Set missing values of soil.df.bulk_ec equal to soil.df.bulk_ec_dc_tc or soil.df.bulk_ec_dc if similar 
    dc_to_non_dc : Calculate missing values of soil.df.bulk_ec based on soil.df.bulk_ec_dc

    """
    if any(np.isnan(soil.df.bulk_ec)):
        BulkECDC(soil)
        conversion(soil)
        dc_to_non_dc(soil)

    return soil.df.bulk_ec.values


def conversion(soil):
    """
    Set missing values of soil.df.bulk_ec equal to soil.df.bulk_ec_dc_tc or soil.df.bulk_ec_dc if similar 

    This function iterates over the soil states to update the 'bulk_ec' values in `soil.df` and
    to annotate corresponding entries in `soil.info['bulk_ec']`. The updates are based on the presence
    of NaN values in the original 'bulk_ec' data and specific conditions related to soil temperature
    and frequency of EC measurements. Annotations indicate the basis for the updated values, referencing
    either `soil.df.bulk_ec_dc_tc` or `soil.df.bulk_ec_dc`, as determined by the conditions met during the
    update process.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: `bulk_ec`, `temperature`, and `frequency_ec`.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - n_states: int
            Number of soil states

    Returns
    -------
    None
        This function updates the soil object in-place and does not return any value.

    Notes
    -----
    - The function checks for NaN values in 'bulk_ec' and updates them based on the soil temperature being
      exactly 298.15K and the frequency of EC measurements being 5Hz or less.
    - Annotations in `soil.info['bulk_ec']` provide insight into the source of the updated values, whether
      they are derived from 'bulk_ec_dc_tc' or 'bulk_ec_dc', aiding in the traceability of the data.
    """
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Equal to soil.df.bulk_ec_dc_tc in predict.bulk_ec.conversion" if np.isnan(soil.df.bulk_ec[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5
                        or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Equal to soil.df.bulk_ec_dc_tc in predict.bulk_ec.conversion"
                        else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec'] = [soil.df.bulk_ec_dc_tc[x] if np.isnan(soil.df.bulk_ec[x]) and soil.df.temperature[x] == 298.15 and soil.df.frequency_ec[x] <= 5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec.conversion" if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] <= 5
                        or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Equal to soil.df.bulk_ec_dc in predict.bulk_ec.conversion"
                        else soil.info.bulk_ec[x] for x in range(soil.n_states)]

    soil.df['bulk_ec'] = [soil.df.bulk_ec_dc[x] if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] <= 5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]


def dc_to_non_dc(soil):
    """
    Calculate missing values of soil.df.bulk_ec based on soil.df.bulk_ec_dc

    This function iterates over the soil states to update the 'bulk_ec' values in `soil.df` for cases where 
    the electromagnetic frequency is 5Hz or higher. It applies the LongmireSmithEC function to account for 
    the frequency shift from zero Hz to the actual frequency. The updated 'bulk_ec' values are then rounded 
    off to a precision defined by `soil.roundn+3`. Corresponding entries in `soil.info['bulk_ec']` are annotated 
    to indicate that the adjustment was made using the LongmireSmithEC function from the `predict.bulk_ec.dc_to_non_dc` 
    module.

    Parameters
    ----------
    soil : object
        A custom soil object containing:

        - df : DataFrame
            Data Frame containing the quantitative information of all soil array-like attributes for each state. 
            Includes: frequency_ec, bulk_ec, bulk_ec_dc and other relevant attributes.
        - info: DataFrame
            Data Frame containing descriptive information about how each array-like attribute was calculated.
        - roundn : int
            Number of decimal places to round results.
        - n_states : int
            Number of soil states.

    Returns
    -------
    None
        This function updates the soil object in-place and does not return any value.

    Notes
    -----
    - The function specifically targets soil states with a frequency of EC measurements of 5Hz or higher.
    - Annotations in `soil.info['bulk_ec']` detail the use of the LongmireSmithEC function for EM frequency 
      adjustments, enhancing the traceability of the data adjustments.

    External Functions
    ------------------
    - LongmireSmithEC : Calculate the soil bulk real electrical conductivity using the Longmire-Smith model and return
    """
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                            if (np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + 
                            "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df["bulk_ec"] = [round(LongmireSmithEC(soil.df.bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) 
                          if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]

