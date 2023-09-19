import numpy as np
from pedophysics import instruments

def FrequencyEC(soil): 
    """
    Return and set missing values of the soil.frequency_ec attribute.

    If any value of the frequency_ec attribute is missing (NaN), 
    it will be set to a default value of 0. Corresponding information 
    in the soil's info dataframe will be updated accordingly.

    Parameters
    ----------
    soil : object
        A custom soil object that contains:

        - frequency_ec : array-like
            Frequency of electric conductivity measurement [Hz]
        - df : DataFrame
            Data Frame containing all the quantitative information of soil array-like attributes for each state
        - info : DataFrame
            Data Frame containing the qualitative information about all array-like soil attributes for each state
        - n_states : int
            Number of states or records in the dataframe.

    Returns
    -------
    np.ndarray
        An array of updated frequency of electric conductivity measurement values

    See Also
    --------
    Inst2FreqC : Function to calculate missing frequency_ec attribute based on soil.instrument

    Example
    -------
    >>> sample = Soil()
    >>> sample.df.frequency_ec
    0   NaN
    Name: frequency_ec, dtype: float64
    >>> FrequencyEC(sample)
    >>> sample.df.frequency_ec
    0    0
    Name: frequency_ec, dtype: float64
    """

    # Check if any value of frequency_ec is missing
    if (np.isnan(soil.df.frequency_ec)).any():
        instruments.Inst2FreqC(soil)

        soil.info['frequency_ec'] = ["Set as 0 Hz (direct current) by default" if np.isnan(soil.df.frequency_ec[x]) or soil.info.frequency_ec[x] == "Set as 0 Hz (direct current) by default"
                                     else soil.info.frequency_ec[x] for x in range(soil.n_states)]
        
        soil.df['frequency_ec'] = [0 if np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    return soil.df.frequency_ec.values