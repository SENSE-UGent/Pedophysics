import numpy as np

def BulkECTC(soil):
    """

    """
    # Check if any value of bulk_perm_inf is missing
    if (np.isnan(soil.df.bulk_ec_tc)).any():

        soil.info['bulk_ec_tc'] = ["Not corrected. Equal to soil.df.bulk_ec" if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] == 298.15 or soil.info.bulk_ec_tc[x] == "Not corrected. Equal to soil.df.bulk_ec"
                                        else soil.info.bulk_ec_tc[x] for x in range(soil.n_states)]

        soil.df['bulk_ec_tc'] = [soil.df.bulk_ec.values[x] if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] == 298.15 else soil.df.bulk_ec_tc[x] for x in range(soil.n_states)] # There is no temperaturre correction

        soil.info['bulk_ec_tc'] = ["Calculated using SheetsHendrickx function" if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] != 298.15 or soil.info.bulk_ec_tc[x] == "Calculated using SheetsHendrickx function"
                                        else soil.info.bulk_ec_tc[x] for x in range(soil.n_states)]
        
        soil.df['bulk_ec_tc'] = [SheetsHendrickx(soil.df.bulk_ec.values[x], soil.df.temperature.values[x]) if np.isnan(soil.df.bulk_ec_tc[x]) and soil.df.temperature[x] != 298.15
                                       else soil.df.bulk_ec_tc[x] for x in range(soil.n_states)]
        
    return soil.df.bulk_ec_tc.values
