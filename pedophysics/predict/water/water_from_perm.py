import numpy as np
from scipy.optimize import minimize

from pedophysics.utils.stats import R2_score

from ..bulk_perm_inf import BulkPermInf
from ..particle_density import ParticleDensity
from ..air_perm import AirPerm
from ..solid_perm import SolidPerm
from ..water_perm import WaterPerm
#from ..cation_exchange_capacity import CEC
from ..texture import Texture

from pedophysics.pedophysical_models.water import LR, LR_W, LR_MV
from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP


############################# Water prediction from permittivity ###################################### 

def WaterFromPerm(soil):
    """ 
    """
    # Condition for constant permittivity frequency
    if np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]):
        fixed_freq(soil)

    # Condition for changing permittivity frequency
    else:
        changing_freq(soil)

    # Converting negative results due to fitting to zero
    soil.info['water'] = [str(soil.info.water[x]) + "--> Set to 0 because of < 0 results" if soil.df.water[x]<0 
                          or soil.info.water[x] ==str(soil.info.water[x]) + "--> Set to 0 because of < 0 results"
                            else soil.info.water[x] for x in range(soil.n_states)]

    soil.df['water'] = [ 0 if soil.df.water[x]<0 else soil.df.water[x] for x in range(soil.n_states)] 


######################################  Changing frequency module #####################################

def changing_freq(soil):    
    '''
    '''
    print('Water perm changing freq')    
    BulkPermInf(soil)    
    bulk_ec = []

    # Defining minimization function to obtain EC
    def objective(bulk_ec, perm_inf, freq, bulk_perm):
        LS_perm = LongmireSmithP(bulk_ec, perm_inf, freq)
        return (LS_perm - bulk_perm)**2

    # Calculating bulk EC from bulk perm when unknown
    for x in range(soil.n_states):
        if np.isnan(soil.df.bulk_ec[x]):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x], soil.bulk_perm[x]), bounds=[(1e-6, 1)])
            bulk_ec.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))
        else:
            bulk_ec.append(np.nan)


    # Saving calculated bulk_ec and its info
    soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.changing_freq" if np.isnan(soil.df.bulk_ec[x])
                            or soil.info.bulk_ec[x] == str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.changing_freq"
                            else soil.info.bulk_ec[x] for x in range(soil.n_states)]
        
    soil.df['bulk_ec'] = [bulk_ec[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]

    # Saving calculated frequency_ec and its info
    soil.info['frequency_ec'] = [str(soil.info.frequency_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.changing_freq" if not np.isnan(bulk_ec[x])
                            or soil.info.frequency_ec[x] == str(soil.info.frequency_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.changing_freq"
                            else soil.info.frequency_ec[x] for x in range(soil.n_states)]
        
    soil.df['frequency_ec'] = [0 if not np.isnan(bulk_ec[x]) else soil.df.frequency_ec[x] for x in range(soil.n_states)]


##################################### Fixed frequnecy ##########################################

def fixed_freq(soil):
    '''
    '''            
    # Condition for fitting approach
    if sum(not np.isnan(soil.water[x]) and not np.isnan(soil.bulk_perm[x]) for x in range(soil.n_states)) >= 3:
        fitting(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_perm[x]) and 5 <= soil.df.frequency_perm[x] <=30e9 for x in range(soil.n_states)):
        non_fitting(soil)


######################################  fixed frequency fitting  #####################################

def fitting(soil):
    '''
    '''
    print('fitting')
    WaterPerm(soil)                   

    # Defining model parameters
    valids = ~np.isnan(soil.df.water) & ~np.isnan(soil.df.bulk_perm) # States where calibration data are
    water_init = np.nanmin(soil.df.water[valids])
    bulk_perm_init = np.nanmin(soil.df.bulk_perm[valids])
    bulk_perm_final = np.nanmax(soil.df.bulk_perm[valids])
    bulk_perm_range = [round(bulk_perm_init - (bulk_perm_final-bulk_perm_init)/soil.range_ratio, soil.roundn), 
                       round(bulk_perm_final + (bulk_perm_final-bulk_perm_init)/soil.range_ratio, soil.roundn)]
    if bulk_perm_range[0] < 0:
        bulk_perm_range[0] = 0
        
    # Obtain Lw attribute if unknown
    if np.isnan(soil.Lw):

        # Defining minimization function to obtain Lw
        def objective_Lw(Lw):
            wund_eval = [WunderlichP(soil.df.water[x], bulk_perm_init, water_init, soil.df.water_perm[x], Lw)[0] if valids[x] else np.nan for x in range(soil.n_states)]
            Lw_RMSE = np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm.values)**2))
            return Lw_RMSE
        
        # Calculating optimal Lw
        result = minimize(objective_Lw, 0.1, bounds=[(-0.2, 0.8)], method='L-BFGS-B')
        soil.Lw = result.x[0]

    # If Lw is known
    if ~np.isnan(soil.Lw):
        if not isinstance(soil.Lw, np.floating):
            soil.Lw = soil.Lw[0]
        Wat_wund = []

        # Defining minimization function to obtain water
        def objective_wat(wat, i):
            return (WunderlichP(wat, bulk_perm_init, water_init, soil.df.water_perm[i], soil.Lw) - soil.df.bulk_perm[i])**2
        
        # Looping over soil states to obtain water using WunderlichP function
        for i in range(soil.n_states):

            if min(bulk_perm_range) <= soil.df.bulk_perm[i] <= max(bulk_perm_range) and ~np.isnan(soil.df.bulk_perm[i]):
                result = minimize(objective_wat, 0.15, args=(i), bounds=[(0, .65)], method='L-BFGS-B')
                Wat_wund.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn))

            else:
                Wat_wund.append(np.nan)

        # Calculating the R2 score of the model fitting
        print("soil.df.water.values", soil.df.water.values)
        print('Wat_wund', Wat_wund)
        R2 = round(R2_score(soil.df.water.values, Wat_wund), soil.roundn)

        # Saving calculated bulk_perm and its info with R2 and valid bulk_perm range
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.water.water_from_perm.fitting, for soil.bulk_perm values between: "+str(bulk_perm_range) if min(bulk_perm_range) <= soil.df.bulk_perm[x] <= max(bulk_perm_range) and ~np.isnan(soil.df.bulk_perm[x]) and np.isnan(soil.df.water[x])
                                or soil.info.water[x] == str(soil.info.water[x]) + "--> Calculated by fitting (R2="+str(R2)+") WunderlichP function in predict.water.water_from_perm.fitting, for soil.bulk_perm values between: "+str(bulk_perm_range)
                                else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [Wat_wund[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]


######################################  fixed frequency non-fitting  #####################################

def non_fitting(soil):
    '''
    '''
    print('non-fitting')
    ParticleDensity(soil)                     
    AirPerm(soil)                      
    SolidPerm(soil)                   
    WaterPerm(soil)              
    Texture(soil)                     
    #CEC(soil)                      

    # Condition for EM frequencies between 5 and 30e6
    if ((soil.df.frequency_perm >= 5) & (soil.df.frequency_perm < 30e6)).all():
        BulkPermInf(soil)

        bulk_ec = []
        # Defining minimization function to obtain bulk_ec using LongmireSmithP
        def objective(bulk_ec, perm_inf, freq_perm, bulk_perm):
            LS_perm = LongmireSmithP(bulk_ec, perm_inf, freq_perm)
            return (LS_perm - bulk_perm)**2
        
        # Calculating bulk_ec
        for i in range(soil.n_states):
            result = minimize(objective, 0.05, args=(soil.df.bulk_perm_inf[i], soil.df.frequency_perm[i], soil.df.bulk_perm[i]), bounds=[(1e-6, 1)], method='L-BFGS-B')
            bulk_ec.append(np.nan if np.isnan(result.fun) else round(result.x[0], soil.roundn+2))

        # Saving calculated bulk_ec and its info
        soil.info['bulk_ec'] = [str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.non_fitting" if np.isnan(soil.df.bulk_ec[x]) 
                                or soil.info.bulk_ec[x] ==str(soil.info.bulk_ec[x]) + "--> Calculated using LongmireSmithP function in predict.water.water_from_perm.non_fitting" else soil.info.bulk_ec[x] for x in range(soil.n_states)]
        
        soil.df['bulk_ec'] = [bulk_ec[x] if np.isnan(soil.df.bulk_ec[x]) else soil.df.bulk_ec[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 30e6 and 100e6
    elif ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR_MV function (reported R2=0.93) in predict.water.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                              or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR_MV function (reported R2=0.93) in predict.water.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]

        soil.df['water'] = [round(LR_MV(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) 
                            if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 100e6 and 200e6
    elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.water.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                        or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR function (reported RMSE=0.032) in predict.water.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
                
        soil.df['water'] = [round(LR(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.alpha), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    # Condition for EM frequencies between 200e6 and 30e9
    elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 30e9))).all():

        # Saving calculated water and its info
        soil.info['water'] = [str(soil.info.water[x]) + "--> Calculated using LR_W function in predict.water.water_from_perm.non_fitting" if np.isnan(soil.df.water[x]) 
                        or soil.info.water[x] ==str(soil.info.water[x]) + "--> Calculated using LR_W function in predict.water.water_from_perm.non_fitting" else soil.info.water[x] for x in range(soil.n_states)]
        
        soil.df['water'] = [round(LR_W(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.air_perm[x], soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)] 
