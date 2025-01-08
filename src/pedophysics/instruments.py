import pandas as pd
import numpy as np

def Inst2FreqP(soil):
    """
    Set missing values of soil.df.frequency_perm and return
    
    This function iterates over each state in the 'soil' object. It checks the instrument type associated with the soil ('GPR', 'TDR', or 'HydraProbe'). 
    If the current 'frequency_perm' value is NaN, the function sets a default frequency value specific to the instrument type 
    (1e9 Hz for 'GPR', 200e6 Hz for 'TDR', and 50e6 Hz for 'HydraProbe') in both the 'info' attribute and the 'df' DataFrame. 
    If the 'frequency_perm' value is already set or the instrument type does not match any specified case, the function retains the existing value.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state

    Returns
    -------
    np.ndarray
        Array containing the updated frequency of dielectric permittivity measurement values

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    """

    soil.info['frequency_perm'] = ["Set as 1e9 Hz because soil.instrument == GPR" if ((soil.instrument == 'GPR') & np.isnan(soil.df.frequency_perm[x])) or (soil.info.frequency_perm[x] == "Set as 1e9 because soil.instrument == GPR")
                                   else soil.info.frequency_perm[x] for x in range(soil.n_states)]

    soil.df['frequency_perm'] = [1e9 if (soil.instrument == 'GPR') & np.isnan(soil.df.frequency_perm[x]) else soil.df.frequency_perm[x] for x in range(soil.n_states)]

    soil.info['frequency_perm'] = ["Set as 200e6 Hz because soil.instrument == TDR" if ((soil.instrument == 'TDR') & np.isnan(soil.df.frequency_perm[x])) or (soil.info.frequency_perm[x] == "Set as 200e6 because soil.instrument == TDR")
                                   else soil.info.frequency_perm[x] for x in range(soil.n_states)]
    
    soil.df['frequency_perm'] = [200e6 if (soil.instrument == 'TDR') & np.isnan(soil.df.frequency_perm[x]) else soil.df.frequency_perm[x] for x in range(soil.n_states)]

    soil.info['frequency_perm'] = ["Set as 50e6 Hz because soil.instrument == HydraProbe" if ((soil.instrument == 'HydraProbe') & np.isnan(soil.df.frequency_perm[x])) or (soil.info.frequency_perm[x] == "Set as 50e6 because soil.instrument == HydraProbe")
                                   else soil.info.frequency_perm[x] for x in range(soil.n_states)]
    
    soil.df['frequency_perm'] = [50e6 if (soil.instrument == 'HydraProbe') & np.isnan(soil.df.frequency_perm[x]) else soil.df.frequency_perm[x] for x in range(soil.n_states)]

    return soil.df.frequency_perm


def Inst2FreqEC(soil):
    """
    Set missing values of soil.df.frequency_ec and return

    This function iterates over each state in the 'soil' object, checking the type of instrument associated with the 'soil' object ('EMI Dualem' or 'EMI EM38-DD'). 
    If the current 'frequency_ec' value is NaN, the function sets a default frequency value (9e3 Hz for 'EMI Dualem' and 16e3 Hz for 'EMI EM38-DD') 
    in both the 'info' attribute and the 'df' dataframe. 
    If the 'frequency_ec' value is already set or the instrument type does not match, the function retains the existing value.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state

    Returns
    -------
    np.ndarray
        Array containing the updated frequency of electric conductivity measurement values

    Notes
    -----
    This function modifies the soil object in-place by updating the `df` and `info` dataframes.
    """
    
    soil.info['frequency_ec'] = ["Set as 9e3 Hz because soil.instrument == EMI Dualem" if (soil.instrument == 'EMI Dualem') & np.isnan(soil.df.frequency_ec[x]) or (soil.info.frequency_ec[x] == "Set as 9e3 because soil.instrument == EMI Dualem")
                                   else soil.info.frequency_ec[x] for x in range(soil.n_states)]
    
    soil.df['frequency_ec'] = [9e3 if (soil.instrument == 'EMI Dualem') & np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    soil.info['frequency_ec'] = ["Set as 16e3 Hz because soil.instrument == EMI EM38-DD" if ((soil.instrument == 'EMI EM38-DD') & np.isnan(soil.df.frequency_ec[x])) or (soil.info.frequency_ec[x] == "Set as 16e3 because soil.instrument == EMI EM38-DD")
                                   else soil.info.frequency_ec[x] for x in range(soil.n_states)]
    
    soil.df['frequency_ec'] = [16e3 if (soil.instrument == 'EMI EM38-DD') & np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    return soil.df.frequency_ec