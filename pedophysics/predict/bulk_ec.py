import numpy as np

from .bulk_ec_dc import BulkECDC

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC

def BulkEC(soil):
    """ 

    """
    if any(np.isnan(soil.df.bulk_ec)):
        conversion(soil)
        BulkECDC(soil)
        dc_to_non_dc(soil)

    return soil.df.bulk_ec.values


def conversion(soil):
    """
    
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
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                            if (np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df["bulk_ec"] = [round(LongmireSmithEC(soil.df.bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) 
                          if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]
    


