from .water_perm import *
from .frequency_perm import *
from .bulk_ec import *
from .bulk_perm_inf import *
from .particle_density import *
from .air_perm import *
from .solid_perm import *

from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP, LR, LR_W, LR_MV
from pedophysics.utils.stats import R2_score

import numpy as np


############################# Bulk permittivity prediction ######################################

def BulkPerm(soil):
    """ 
    """
    if (np.isnan(soil.df.bulk_perm)).any():  # Go over if any value is missing        
        FrequencyPerm(soil)

        # Condition to ask for frequency data
        if (np.isnan(soil.df.frequency_perm)).all():
            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Unmodified value. Please provide soil.frequency_perm" if True else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            soil.df['bulk_perm'] = [soil.df.bulk_perm[x] if True else soil.df.bulk_perm[x] for x in range(soil.n_states)]

        # Condition for fixed EM frequency
        elif np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]):
            fixed_freq(soil)

        # Condition for changing EM frequency
        else:
            changing_freq(soil)

    return soil.df.bulk_perm.values


############################################ Fixed frequency ##########################################

def fixed_freq(soil):
    '''
    
    '''
    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_perm[x]) for x in range(soil.n_states)) >= 3 :
        fitting(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.bulk_perm[x]) and not np.isnan(soil.df.water[x]) and soil.df.frequency_perm[x] >= 5 and soil.df.frequency_perm[x] <= 30e9 for x in range(soil.n_states) ):
        non_fitting(soil)

        
########################################### Fitting routine ##################################################

def fitting(soil):
    '''
    '''
    Temperature(soil)
    WaterPerm(soil)                      

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_perm) # States where calibration data are
    water_init = min(soil.df.water[valids])
    bulk_perm_init = min(soil.df.bulk_perm[valids])
    water_final = max(soil.df.water[valids])
    water_range = [round(water_init - (water_final-water_init)/soil.range_ratio, soil.roundn), 
                  round(water_final + (water_final-water_init)/soil.range_ratio, soil.roundn)]
    if water_range[0] < 0:
        water_range[0] = 0

    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]    
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm)**2))
            return Lw_RMSE
    
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        # Calculating the R2 score of the model fitting
        R2 = round(R2_score(soil.df.bulk_perm, WunderlichP(soil.df.water, bulk_perm_init, water_init, soil.df.water_perm, soil.Lw)), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid water range
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.bulk_perm.fitting, for soil.water values between"+str(water_range) if ((min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_perm[x]))
                                or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.bulk_perm.fitting, for soil.water values between"+str(water_range)
                                else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
        soil.df['bulk_perm'] = [round(WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], soil.Lw), soil.roundn) 
                              if (min(water_range) <= soil.water[x] <= max(water_range)) and np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]


############################################# Non fitting ###############################################

def non_fitting(soil):
    '''
    
    '''
    # Condition for lowest EM frequency
    if any(np.isnan(soil.df.bulk_perm[x]) and soil.df.frequency_perm[x] >= 5 and soil.df.frequency_perm[x] < 30e6 for x in range(soil.n_states)):

        BulkPermInf(soil)              
        BulkEC(soil)

        # Saving calculated bulk_perm and its info
        soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.non_fitting"
                                else soil.info.bulk_perm[x] for x in range(soil.n_states)]
        
        soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_ec[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) 
                                if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

    # Condition for EM frequency of common moisture sensors and GPR
    elif (np.isnan(soil.df.bulk_perm)).any() & ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm <= 30e9)).all(): 
        Temperature(soil)
        ParticleDensity(soil)                      
        AirPerm(soil)                   
        SolidPerm(soil)                  
        WaterPerm(soil)               
        Texture(soil)                    

        if ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():

            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR_MV (reported R2=0.93) function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR_MV (reported R2=0.93) function in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR_MV(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

        elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

            if np.isnan(soil.alpha): soil.alpha = 0.5 
            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

        elif ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 30e9)).all(): 

            soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LR_W function in predict.bulk_perm.non_fitting" if np.isnan(soil.df.bulk_perm[x])
                                    or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LR_W function in predict.bulk_perm.non_fitting"
                                    else soil.info.bulk_perm[x] for x in range(soil.n_states)]
            
            soil.df['bulk_perm'] = [np.round(LR_W(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) 
                                    if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)] 


########################################### Changing frequency ############################################

def changing_freq(soil):
    '''
    '''
    BulkPermInf(soil)             
    BulkEC(soil)

    # Saving calculated bulk_perm and its info
    soil.info['bulk_perm'] = [str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.changing_freq" if np.isnan(soil.df.bulk_perm[x])
                            or soil.info.bulk_perm[x] == str(soil.info.bulk_perm[x]) + "--> Calculated using LongmireSmithP function in predict.bulk_perm.changing_freq"
                            else soil.info.bulk_perm[x] for x in range(soil.n_states)]
    
    soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_ec[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) 
                            if (np.isnan(soil.df.bulk_perm[x])) else soil.df.bulk_perm[x] for x in range(soil.n_states)]