import numpy as np

def Inst2FreqP(soil):
    """ 

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

    """
    soil.info['frequency_ec'] = ["Set as 9e3 Hz because soil.instrument == EMI Dualem" if (soil.instrument == 'EMI Dualem') & np.isnan(soil.df.frequency_ec[x]) or (soil.info.frequency_ec[x] == "Set as 9e3 because soil.instrument == EMI Dualem")
                                   else soil.info.frequency_ec[x] for x in range(soil.n_states)]
    
    soil.df['frequency_ec'] = [9e3 if (soil.instrument == 'EMI Dualem') & np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    soil.info['frequency_ec'] = ["Set as 16e3 Hz because soil.instrument == EMI EM38-DD" if ((soil.instrument == 'EMI EM38-DD') & np.isnan(soil.df.frequency_ec[x])) or (soil.info.frequency_ec[x] == "Set as 16e3 because soil.instrument == EMI EM38-DD")
                                   else soil.info.frequency_ec[x] for x in range(soil.n_states)]
    
    soil.df['frequency_ec'] = [16e3 if (soil.instrument == 'EMI EM38-DD') & np.isnan(soil.df.frequency_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]

    return soil.df.frequency_ec