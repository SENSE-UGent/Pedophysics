import numpy as np

def Inst2FreqP(soil):
    """
    Return or compute missing values of the soil.frequency_perm attribute.

    If the `frequency_perm` value of the soil is missing (`np.nan`), this function assigns 
    a default frequency based on the type of instrument used to measure the soil:
    - GPR: 1e9 Hz
    - TDR: 200e6 Hz
    - HydraProbe: 50e6 Hz

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - frequency_perm : array-like
            Frequency of dielectric permittivity measurement [Hz]
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


def Inst2FreqC(soil):
    """
    Return or compute missing values of the soil.frequency_ec attribute.

    If the `frequency_ec` value of the soil is missing (`np.nan`), this function assigns 
    a default frequency based on the type of instrument used to measure the soil:
    - EMI Dualem: 9e3 Hz
    - EMI EM38-DD: 16e3 Hz

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - frequency_ec : array-like
            Frequency of electric conductivity measurement [Hz]
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