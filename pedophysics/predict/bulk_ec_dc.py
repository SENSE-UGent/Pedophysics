import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC, SheetsHendrickx


def BulkECDC(soil):
    """
    Compute missing values of soil.df.bulk_ec_dc and return

    This function checks for NaN values in `soil.df.bulk_ec_dc`. If NaN values are found, it attempts to estimate
    these missing values using a series of steps. First, it calls the `BulkECDCTC` function to potentially fill in
    missing `bulk_ec_dc` values. If after this step there are still NaN values in `bulk_ec_dc` but corresponding
    `bulk_ec_dc_tc` values are available, it calls `tc_to_non_tc` to convert temperature-corrected values to their
    non-temperature-corrected counterparts. Similarly, if there are NaN values in `bulk_ec_dc` but non-direct current
    (non-dc) `bulk_ec` values are available, it calls `non_dc_to_dc` to convert these to direct current values.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            Data Frame containing quantitative information of soil attributes for each state, including 
            `bulk_ec_dc`, `bulk_ec_dc_tc`, and `bulk_ec`.
        - n_states: int
            Number of soil states.

    Returns
    -------
    np.ndarray
        soil.df.bulk_ec_dc.values: soil bulk real electrical conductivity direct current values.

    External functions
    --------
    BulkECDCTC : Compute missing values of soil.df.bulk_ec_dc_tc and return
    tc_to_non_tc : Calculate missing values of soil.df.bulk_ec_dc based on soil.df.bulk_ec_dc_tc
    non_dc_to_dc : Calculate missing values of soil.df.bulk_ec_dc based on soil.df.bulk_ec

    Notes
    -----
    - The function operates in-place, modifying the soil object's DataFrame directly.
    - It relies on the availability of either `bulk_ec_dc_tc` or `bulk_ec` values to estimate missing `bulk_ec_dc` values.
    - The process involves a sequence of conversions and estimations, which may depend on external functions not
      detailed here (`BulkECDCTC`, `tc_to_non_tc`, and `non_dc_to_dc`).
    """    
    from .bulk_ec_dc_tc import BulkECDCTC # Lazy import to avoid circular dependency
    if any(np.isnan(soil.df.bulk_ec_dc)):
        BulkECDCTC(soil)

        if any(np.isnan(soil.df.bulk_ec_dc[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) for x in range(soil.n_states)):
            tc_to_non_tc(soil)

        if any(np.isnan(soil.df.bulk_ec_dc[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)):
            non_dc_to_dc(soil)

    return soil.df.bulk_ec_dc.values


def tc_to_non_tc(soil):
    """
    Calculate missing values of soil.df.bulk_ec_dc based on soil.df.bulk_ec_dc_tc

    This function iterates through soil states to update `bulk_ec_dc` in `soil.df` where it is NaN. For soil states
    at standard temperature (298.15K), `bulk_ec_dc` is set directly equal to `bulk_ec_dc_tc`. For other temperatures,
    it uses the `SheetsHendrickx` function within a minimization process to estimate `bulk_ec_dc` from `bulk_ec_dc_tc`.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A DataFrame containing quantitative information of soil attributes for each state, including `bulk_ec_dc`,
            `bulk_ec_dc_tc`, and `temperature`.
        - info: DataFrame
            A DataFrame containing descriptive information for each soil state, including `bulk_ec_dc`.
        - n_states: int
            The number of soil states.
        - roundn: int
            Number of decimal places to round results.

    Returns
    -------
    None
        This function updates the soil object in-place and does not return any value.

    External functions
    --------
    SheetsHendrickx : Calculate the soil bulk real electrical conductivity using the Sheets-Hendricks model and return

    Notes
    -----
    - The function uses a minimization process to estimate `bulk_ec_dc` from `bulk_ec_dc_tc` for soil states
      not at standard temperature. The objective function aims to minimize the difference between the
      temperature-corrected value obtained using `SheetsHendrickx` and the given `bulk_ec_dc_tc`.
    - It directly sets `bulk_ec_dc` equal to `bulk_ec_dc_tc` for soil states at standard temperature (298.15K)
      without any correction.
    - Updates and calculations are logged in `soil.info` for traceability.
    """    

    # Defining minimization function to obtain DC bulk EC 
    def objective_tc_to_non_tc(bulk_ec_dc, bulk_ec_dc_tc, temperature):
        return (SheetsHendrickx(bulk_ec_dc, temperature) - bulk_ec_dc_tc)**2

    for i in range(soil.n_states):
        if soil.df.temperature[i] == 298.15 and np.isnan(soil.df.bulk_ec_dc[i]):
            soil.info.loc[i, 'bulk_ec_dc'] = str(soil.info.bulk_ec_dc[i]) + "--> Equal to soil.df.bulk_ec_dc_tc in predict.bulk_ec_dc.tc_to_non_tc"
            soil.df.loc[i, 'bulk_ec_dc'] = soil.df.bulk_ec_dc_tc[i]

        elif soil.df.temperature[i] != 298.15 and np.isnan(soil.df.bulk_ec_dc[i]):
            res = minimize(objective_tc_to_non_tc, 0.05, args=(soil.df.bulk_ec_dc_tc[i], soil.df.temperature[i]), bounds=[(0, 1)])

            soil.info.loc[i, 'bulk_ec_dc'] = str(soil.info.bulk_ec_dc[i]) + "--> Calculated from soil.df.bulk_ec_dc_tc using SheetsHendrickx function in predict.bulk_ec_dc.tc_to_non_tc"
            soil.df.loc[i, 'bulk_ec_dc'] = np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2)


def non_dc_to_dc(soil):
    """
    Calculate missing values of soil.df.bulk_ec_dc based on soil.df.bulk_ec

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
        - info : DataFrame
            Data Frame containing descriptive information about how each array-like attribute was determined or modified.

    Notes
    -----
    The function modifies the soil object in-place by updating the `info` attribute to include details about the 
    conversion method used for each state.

    External Functions
    ------------------
    - LongmireSmithEC : Calculate the soil bulk real electrical conductivity using the Longmire-Smith model and return
    """

    # Defining minimization function to obtain DC bulk EC 
    def objective_non_dc_to_dc(x, frequency_ec, bulk_ec):
        return (LongmireSmithEC(x, frequency_ec) - bulk_ec)**2

    for i in range(soil.n_states):
        if (soil.df.frequency_ec[i] <= 5) and np.isnan(soil.df.bulk_ec_dc[i]):
            soil.info.loc[i, 'bulk_ec_dc'] = str(soil.info.bulk_ec_dc[i]) + "--> Equal to soil.df.bulk_ec in predict.bulk_ec_dc.non_dc_to_dc" 
            soil.df.loc[i, 'bulk_ec_dc'] = soil.df.bulk_ec[i]

        elif soil.df.frequency_ec[i] > 5 and np.isnan(soil.df.bulk_ec_dc[i]):
            res = minimize(objective_non_dc_to_dc, 0.05, args=(soil.df.frequency_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 1)])

            soil.info.loc[i, 'bulk_ec_dc'] = str(soil.info.bulk_ec_dc[i]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec_dc.non_dc_to_dc"
            soil.df.loc[i, 'bulk_ec_dc'] = np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2)

