import numpy as np
from scipy.optimize import minimize

from .water_ec import *
from .frequency_ec import *
from .particle_density import *
from .solid_ec import *
from .bulk_ec_dc_tc import *

from pedophysics.pedophysical_models.bulk_ec import WunderlichEC, LongmireSmithEC, Fu, SheetsHendrickx


def BulkECDC(soil):
    """
    
    """        
    BulkECDCTC(soil)
    tc_to_non_tc(soil)

    return soil.df.bulk_ec_dc.values

#    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> Assumed equal to soil.df.bulk_ec because soil.df.frequency_ec[x] >= 5" if any(np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) and soil.df.frequency_ec[x] >= 5) or 
#            soil.info.bulk_ec_dc[x] == str(soil.info.bulk_ec_dc[x]) + "--> Assumed equal to soil.df.bulk_ec because soil.df.frequency_ec[x] >= 5" else soil.info.bulk_ec_dc[x] for x in range(soil.n_states)]

#    soil.df['bulk_ec_dc'] = soil.df.bulk_ec


def tc_to_non_tc(soil):
    """

    """    
    # Defining minimization function to obtain DC bulk EC 
    def objective_func_ec_tc(bulk_ec_dc, bulk_ec_dc_tc, temperature):
        return (SheetsHendrickx(bulk_ec_dc, temperature) - bulk_ec_dc_tc)**2
    bulk_ec_dc = []

    for i in range(soil.n_states):
        res = minimize(objective_func_ec_tc, 0.05, args=(soil.df.bulk_ec_dc_tc[i], soil.df.frequency_ec[i]), bounds=[(0, 1)])
        bulk_ec_dc.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2) )

    soil.info['bulk_ec_dc'] = [str(soil.info.bulk_ec_dc[x]) + "--> Calculated from soil.df.bulk_ec_dc_tc using SheetsHendrickx function in predict.bulk_ec.tc_to_non_tc" if 
                    np.isnan(soil.df.bulk_ec_dc[x]) and not np.isnan(soil.df.bulk_ec_dc_tc[x]) or 
                    soil.info.bulk_ec_dc[x] == str(soil.info.bulk_ec_dc[x]) + "--> Calculated from soil.df.bulk_ec_dc_tc using SheetsHendrickx function in predict.bulk_ec.tc_to_non_tc"  
                    else soil.info.bulk_ec_dc[x] for x in range(soil.n_states)]
    
    soil.df['bulk_ec_dc'] = [round(bulk_ec_dc[i], soil.roundn+3) if np.isnan(soil.df.bulk_ec_dc[i]) and not np.isnan(soil.df.bulk_ec_dc_tc[i]) else soil.df.bulk_ec_dc[i] for i in range(soil.n_states)]


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


