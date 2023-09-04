import numpy as np
from scipy.optimize import minimize

from pedophysics.utils.stats import R2_score

from pedophysics.predict.water_ec import WaterEC
from pedophysics.predict.frequency_ec import *
from pedophysics.predict.particle_density import *
from pedophysics.predict.solid_ec import *
#from pedophysics.predict.cation_exchange_capacity import *

from pedophysics.pedophysical_models.bulk_ec import LongmireSmithEC, Fu, WunderlichEC


############################# Water prediction from electrical ecuctivity ###################################### 

def WaterFromEC(soil):
    """ 
        Soil volumetric water content prediction from real bulk electrical ecuctivity
    """
    FrequencyEC(soil)
    # Condition for non-DC frequency
    print('WaterFromEC')
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 for x in range(soil.n_states)):
        bulk_ec_dc = non_dc_to_dc(soil) 

    # Condition for DC frequency
    else:
        bulk_ec_dc = soil.df.bulk_ec

    print('bulk_ec_dc', bulk_ec_dc)
    dc_freq(soil, bulk_ec_dc)


######################################    non DC frequency   #####################################

# This is a special function to be used in non_dc_freq function
def non_dc_to_dc(soil):
    '''
    '''
    # Defining minimization function to obtain DC bulk EC 
    def objective_func_ec_dc(bulk_ec_dc, frequency_ec, bulk_ec):
        return (LongmireSmithEC(bulk_ec_dc, frequency_ec) - bulk_ec)**2
    ec_dc = []

    for i in range(soil.n_states):
        if soil.df.frequency_ec[i] <= 5 or np.isnan(soil.df.bulk_ec[i]):
            ec_dc.append(soil.df.bulk_ec[i])

        elif soil.df.frequency_ec[i] > 5 and not np.isnan(soil.df.bulk_ec[i]):
            res = minimize(objective_func_ec_dc, 0.05, args=(soil.df.frequency_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 1)])
            ec_dc.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn+2) )

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift done using LongmireSmithEC function in predict.water.water_from_ec.non_dc_freq_to_dc_freq" if 
                    soil.df.frequency_ec[x] > 5 and not np.isnan(soil.df.bulk_ec[x]) or 
                    soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift done using LongmireSmithEC function in predict.water.water_from_ec.non_dc_freq_to_dc_freq" 
                    else soil.info.bulk_ec[x] for x in range(soil.n_states)]
            
    ec_dc = np.array(ec_dc)
    return ec_dc


######################################    DC frequency   #####################################

def dc_freq(soil, bulk_ec_dc):
    '''
    '''
    print('dc_freq')
    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states)) >= 3:
        fitting(soil, bulk_ec_dc)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states)):
        non_fitting(soil, bulk_ec_dc)


######################################  DC frequency - non fitting  #####################################

def non_fitting(soil, bulk_ec_dc):
    '''
    '''
    print('non_fitting')
    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water using Fu
    def objective_func_wat(x, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec, EC):
        return (Fu(x, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec) - EC)**2
    wat = []

    # Calculating water
    for i in range(soil.n_states):
        res = minimize(objective_func_wat, 0.15, args=(soil.df.clay[i], soil.df.bulk_density[i], soil.df.particle_density[i], soil.df.water_ec[i], soil.df.solid_ec[i], 
                                                        soil.df.dry_ec[i], soil.df.sat_ec[i], bulk_ec_dc[i]), bounds=[(0, .65)] )
        wat.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

   # Saving calculated water and its info
    soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water.water_from_ec.non_fitting" if  np.isnan(soil.df.water[x]) or 
                          soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water.water_from_ec.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
    
    soil.df['water'] = [round(wat[i], soil.roundn) if np.isnan(soil.df.water[i]) else soil.df.water[i] for i in range(soil.n_states) ]


######################################  DC fitting  #####################################

def fitting(soil, bulk_ec_dc):
    '''
    '''
    WaterEC(soil) 
    print('water EC fitting')
    print('bulk_ec_dc', bulk_ec_dc)
    print('soil.df.frequency_ec', soil.df.frequency_ec)
    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(bulk_ec_dc) # States where calibration data are
    water_init = np.nanmin(soil.df.water[valids])
    bulk_ec_init = np.nanmin(bulk_ec_dc[valids])
    bulk_ec_final = np.nanmax(bulk_ec_dc[valids])
    bulk_ec_range = [round(bulk_ec_init - (bulk_ec_final-bulk_ec_init)/soil.range_ratio, soil.roundn), 
                     round(bulk_ec_final + (bulk_ec_final-bulk_ec_init)/soil.range_ratio, soil.roundn)]
    if bulk_ec_range[0] < 0:
        bulk_ec_range[0] = 0

    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain water
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - bulk_ec_dc)**2))
            return Lw_RMSE

        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]
        
    # If Lw is known
    if ~np.isnan(soil.Lw):
        Wat_wund = []

        # Defining minimization function to obtain water
        def objective_wat(wat, i):
            Wat_RMSE = np.sqrt((WunderlichEC(wat, bulk_ec_init, water_init, soil.df.water_ec[i], soil.Lw) - bulk_ec_dc[i])**2)
            return Wat_RMSE
        
        # Looping over soil states to obtain water using WunderlichEC function
        for i in range(soil.n_states):
            if (min(bulk_ec_range) <= bulk_ec_dc[i] <= max(bulk_ec_range)) & ~np.isnan(bulk_ec_dc[i]):
                result = minimize(objective_wat, 0.15, args=(i), bounds=[(0, .65)], method='L-BFGS-B')
                Wat_wund.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn))

            else:
                Wat_wund.append(np.nan)

        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.water[valids], np.array(Wat_wund)[valids]), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid bulk_ec range
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.water.water_from_ec.fitting, for soil.bulk_ec values between: "+str(bulk_ec_range) 
                              if min(bulk_ec_range) <= bulk_ec_dc[x] <= max(bulk_ec_range) and np.isnan(soil.df.water[x])
                                or soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.water.water_from_ec.fitting, for soil.bulk_ec values between: "+str(bulk_ec_range)
                                else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [Wat_wund[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Converting negative results due to fitting to zero
    soil.info['water'] = [str(soil.info.water[x]) + "--> Set to 0 because of < 0 results" if soil.df.water[x]<0 or soil.info.water[x] == str(soil.info.water[x]) + "--> Set to 0 because of < 0 results"
                            else soil.info.water[x] for x in range(soil.n_states)]
    
    soil.df['water'] = [ 0 if soil.df.water[x]<0 else soil.df.water[x] for x in range(soil.n_states)] 
