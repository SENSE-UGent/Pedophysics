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

        if any(soil.df.frequency_ec[x] < 5 and np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) for x in range(soil.n_states)):
            dc_freq(soil)

        if any(soil.df.frequency_ec[x] >= 5 and np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) for x in range(soil.n_states)):
            non_dc_freq(soil)

    return soil.df.bulk_ec.values


##########################################    non DC frequency   ###########################################

def non_dc_freq(soil):
    '''
    
    '''
    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    bulk_ec_dc = [round(Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]), soil.roundn+3) 
                  if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]

#    frequency_ec = [0 if np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x] >= 5 else soil.df.frequency_ec[x] for x in range(soil.n_states)]    

    # Calculating bulk EC non-DC using LongmireSmithEC function and save the information
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithEC function in predict.bulk_ec.non_dc_freq" if (soil.df.frequency_ec[x]>=5 and np.isnan(soil.df.bulk_ec[x])) or 
                            soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithEC function in predict.bulk_ec.non_dc_freq"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]

#    soil.df['bulk_ec'] = [round(LongmireSmithEC(bulk_ec_dc[x], soil.df.frequency_ec[x]), soil.roundn+3) if (frequency_ec[x] != soil.df.frequency_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]


##########################################    DC frequency   ###########################################

def dc_freq(soil):
    '''
    
    '''
    # Condition for fitting routine 
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_ec[x]) and soil.df.frequency_ec[x]<5 for x in range(soil.n_states))>= 3:
        fitting(soil)

    # Condition for non-fitting routine 
    if any(not np.isnan(soil.df.water[x]) and np.isnan(soil.df.bulk_ec[x]) and soil.df.frequency_ec[x]<5 for x in range(soil.n_states)):
        non_fitting(soil)


######################################  DC frequency - fitting  ##########################################

def fitting(soil):
    '''
    
    '''
    WaterEC(soil)                    

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_ec) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_ec_init = min(soil.df.bulk_ec[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], Lw)[0] if valids[x] and soil.df.frequency_ec[x]<5 else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_ec)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):

        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.bulk_ec, WunderlichEC(soil.df.water, bulk_ec_init, water_init, soil.df.water_ec, soil.Lw)), soil.roundn)

        # Saving calculated bulk_ec and its info with R2 and valid water range
        soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_ec[x]))
                                or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichEC function in predict.bulk_ec.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_ec[x] for x in range(soil.n_states)]
                
        soil.df['bulk_ec'] = [round(WunderlichEC(soil.df.water[x], bulk_ec_init, water_init, soil.df.water_ec[x], soil.Lw), soil.roundn+3) 
                              if (min(water_range) <= soil.water[x] <= max(water_range)) and soil.df.frequency_ec[x]<5 and np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]


######################################  DC frequency - non fitting  #####################################

def non_fitting(soil):
    '''
    
    '''
    Texture(soil)
    ParticleDensity(soil)
    WaterEC(soil)
    SolidEC(soil)

    # Calculating bulk EC DC using Fu function and save the information
    bulk_ec_dc = [Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_ec[x], soil.df.solid_ec[x], soil.df.dry_ec[x], soil.df.sat_ec[x]) 
                  if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]   

    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting" if (soil.df.frequency_ec[x]<5 and np.isnan(soil.df.bulk_ec[x])) 
                            or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.bulk_ec.non_fitting"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
 
    soil.df['bulk_ec'] = [round(bulk_ec_dc[x], soil.roundn+3) if soil.df.frequency_ec[x]<5 else soil.df.bulk_ec[x] for x in range(soil.n_states)]
