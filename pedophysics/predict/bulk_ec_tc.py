import numpy as np

from .bulk_ec import BulkEC

from pedophysics.pedophysical_models.bulk_ec import SheetsHendrickx

def BulkECTC(soil):
    """
    Calculate missing values of soil.df.bulk_ec_tc and return

    This function checks if any values of `bulk_ec_tc` (temperature-corrected bulk electrical conductivity)
    in `soil.df` are missing. If so, it updates these values based on the soil's temperature. For temperatures
    exactly at 298.15K, `bulk_ec_tc` is set equal to `bulk_ec` without correction. For other temperatures,
    `bulk_ec_tc` is calculated using the SheetsHendrickx function. Corresponding entries in `soil.info['bulk_ec_tc']`
    are annotated to indicate whether values were directly taken from `bulk_ec` or calculated using the
    SheetsHendrickx function.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing quantitative information of soil attributes for each state, including 
            `bulk_ec`, `bulk_ec_tc`, and `temperature`.
        - info: DataFrame
            Dictionary containing descriptive information about how each attribute was calculated.
        - n_states: int
            Number of soil states.

    Returns
    -------
    np.ndarray
        soil.df.bulk_ec_tc.values: an array of updated Soil bulk real electrical conductivity temperature corrected (298.15 K) values

    External functions
    --------
    BulkEC : Calculate missing values of soil.df.bulk_ec and return
    SheetsHendrickx : Calculate the soil bulk real electrical conductivity using the Sheets-Hendricks model and return

    Notes
    -----
    - If `bulk_ec_tc` is NaN and temperature is 298.15K, `bulk_ec_tc` is set equal to `bulk_ec` with an annotation
      indicating no temperature correction was applied.
    - If `bulk_ec_tc` is NaN and temperature is not 298.15K, `bulk_ec_tc` is calculated using the SheetsHendrickx
      function, with an annotation indicating it was temperature-corrected.
    - The function updates the soil object in-place and returns an array of the updated `bulk_ec_tc` values.

    """
    # Check if any value of bulk_perm_inf is missing
    if (np.isnan(soil.df.bulk_ec_tc)).any():
        BulkEC(soil)

        soil.info['bulk_ec_tc'] = [str(soil.info.bulk_ec_tc[x]) + "--> Not corrected. Equal to soil.df.bulk_ec" if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] == 298.15 or 
                                   soil.info.bulk_ec_tc[x] == str(soil.info.bulk_ec_tc[x]) + "--> Not corrected. Equal to soil.df.bulk_ec" else soil.info.bulk_ec_tc[x] for x in range(soil.n_states)]

        soil.df['bulk_ec_tc'] = [soil.df.bulk_ec.values[x] if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] == 298.15 else soil.df.bulk_ec_tc[x] for x in range(soil.n_states)] # There is no temperaturre correction

        soil.info['bulk_ec_tc'] = [str(soil.info.bulk_ec_tc[x]) + "--> Calculated using SheetsHendrickx function" if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] != 298.15 or 
                                   soil.info.bulk_ec_tc[x] == str(soil.info.bulk_ec_tc[x]) + "--> Calculated using SheetsHendrickx function" else soil.info.bulk_ec_tc[x] for x in range(soil.n_states)]
        
        soil.df['bulk_ec_tc'] = [SheetsHendrickx(soil.df.bulk_ec.values[x], soil.df.temperature.values[x]) if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] != 298.15
                                       else soil.df.bulk_ec_tc[x] for x in range(soil.n_states)]
        
    return soil.df.bulk_ec_tc.values
