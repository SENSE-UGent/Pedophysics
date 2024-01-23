import numpy as np
from scipy.optimize import minimize

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC, SheetsHendrickx


def BulkECDC(soil):
    """
    
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

