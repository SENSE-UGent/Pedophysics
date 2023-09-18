import numpy as np
from scipy.optimize import minimize

from .water_ec import *
from .frequency_ec import *
from .particle_density import *
from .solid_ec import *

from pedophysics.pedophysical_models.bulk_ec import WunderlichEC, LongmireSmithEC, Fu


######################################  Predict Bulk Electrical Conductivity   #####################################

def BulkEC(soil):
    '''

    '''      
    if (np.isnan(soil.df.bulk_ec)).any():  # Go over if any value is missing        
        FrequencyEC(soil)

        if any(soil.df.frequency_ec[x] >= 5 and np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) for x in range(soil.n_states)):
            bulk_ec_dc = non_dc_to_dc(soil)

        else:
            bulk_ec_dc = soil.df.bulk_ec

        bulk_ec_dc = dc_freq(soil, bulk_ec_dc)
        dc_to_non_dc(soil, bulk_ec_dc)

    return soil.df.bulk_ec.values

##########################################    non DC to DC frequency   ###########################################


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

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" if 
                    soil.df.frequency_ec[x] > 5 and not np.isnan(soil.df.bulk_ec[x]) or 
                    soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from actual to zero Hz using LongmireSmithEC function in predict.bulk_ec.non_dc_to_dc" 
                    else soil.info.bulk_ec[x] for x in range(soil.n_states)]
            
    return np.array(ec_dc)

##########################################    DC frequency   ###########################################


def dc_freq(soil, bulk_ec_dc):
    '''
    
    '''
    # Condition for fitting routine 
    if sum(not np.isnan(soil.water[x]) and not np.isnan(bulk_ec_dc[x]) for x in range(soil.n_states))>= 3:
        bulk_ec_dc = fitting(soil, bulk_ec_dc)

    # Condition for non-fitting routine 
    if any(not np.isnan(soil.df.water[x]) and np.isnan(bulk_ec_dc[x])  for x in range(soil.n_states)):
        bulk_ec_dc = non_fitting(soil, bulk_ec_dc)

    return bulk_ec_dc


######################################  DC frequency - fitting  ##########################################

def fitting(soil, bulk_ec_dc):
    '''
    
    '''
    WaterEC(soil)                    

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(bulk_ec_dc) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_ec_init = min(bulk_ec_dc[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - bulk_ec_dc)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(bulk_ec_dc, WunderlichEC(soil.df.water, bulk_ec_init, water_init, soil.df.water_ec, soil.Lw)), soil.roundn)

        # Saving calculated bulk_ec and its info with R2 and valid water range
        soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_ec[x]))
                                or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
                
        bulk_ec_dc = [round(WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], soil.Lw), soil.roundn+3) if 
                      (min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(bulk_ec_dc[x]) else bulk_ec_dc[x] for x in range(soil.n_states)]
       
        return bulk_ec_dc


######################################  DC frequency - non fitting  #####################################

def non_fitting(soil, bulk_ec_dc):
    '''
    
    '''
    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    # Calculating bulk EC DC using Fu function and save the information
    #bulk_ec_dc = [Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]) 
    #              if np.isnan(bulk_ec_dc[x]) else bulk_ec_dc[x] for x in range(soil.n_states)]   

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting" if np.isnan(bulk_ec_dc[x]) 
                            or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
 
    bulk_ec_dc = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) if np.isnan(bulk_ec_dc[x]) else bulk_ec_dc[x] for x in range(soil.n_states)]

    return bulk_ec_dc


##########################################  DC to non DC frequency   ###########################################

def dc_to_non_dc(soil, bulk_ec_dc):
    '''

    '''   
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" if (np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or 
                soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> EM frequency shift from zero Hz to actual using LongmireSmithEC function in predict.bulk_ec.dc_to_non_dc" 
                else soil.info.bulk_ec[x] for x in range(soil.n_states)]

    bulk_ec_non_dc = [round(LongmireSmithEC(bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else bulk_ec_dc[x] for x in range(soil.n_states)]

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" if (not np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5) or 
                soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Set to value given by the user" 
                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
    
    soil.df["bulk_ec"] = [bulk_ec_non_dc[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]