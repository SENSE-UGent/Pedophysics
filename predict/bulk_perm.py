import numpy as np
from predict.frequency_perm import frequency_perm
from predict.water_perm import water_perm
from predict.bulk_perm_inf import bulk_perm_inf

from predict.particle_density import particle_density            
from predict.air_perm import air_perm      
from predict.solid_perm import solid_perm 
from predict.CEC import CEC
from predict.texture import texture

import pedophysical_models as pm


def bulk_perm(soil):
    """ 
        Real bulk dielectric permittivity prediction
    """
    if (np.isnan(soil.bulk_perm) == True).any:  # Go over if any value is missing


        ################################################### Fitting routine ##################################################
        soil.frequency_perm = frequency_perm(soil)

        if (np.all(soil.frequency_perm == soil.frequency_perm[0])): # No changes in frequency
            
            if ((sum([1 if ( (np.isnan(soil.water[x]) == False) & (np.isnan(soil.bulk_perm[x]) == False) ) else 0 for x in range(soil.n_states)]) >= 3)) :
                soil.water_perm = water_perm(soil)                        # Yellow

                if (np.isnan(soil.Lw == True) or np.isnan(soil.water_init == True) or np.isnan(soil.bulk_perm_init == True)).any:
                    soil.water_init = min([soil.df.water[x] if ( (np.isnan(soil.df.water[x]) == False) & (np.isnan(soil.df.bulk_perm[x]) == False) ) else np.nan for x in range(soil.n_states)])  
                    soil.bulk_perm_init = min([soil.df.bulk_perm[x] if ( (np.isnan(soil.df.water[x]) == False) & (np.isnan(soil.df.bulk_perm[x]) == False) ) else np.nan for x in range(soil.n_states)])
                    water_final = max([soil.df.water[x] if ( (np.isnan(soil.df.water[x]) == False) & (np.isnan(soil.df.bulk_perm[x]) == False) ) else np.nan for x in range(soil.n_states)])
                    water_range = [soil.water_init - (water_final-soil.water_init)/5, water_final + (water_final-soil.water_init)/5]

                    Lw_ = np.arange(-0.2,  1,  (1.2) / 100 )
                    Wund_RMSE = []

                    for i in range(len(Lw_)):
                        wund_eval = [pm.wunderlich_p(soil.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], Lw_[i]) if ( (np.isnan(soil.water[x]) == False) & (np.isnan(soil.bulk_perm[x]) == False) ) else np.nan for x in range(soil.n_states)]
                        Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.bulk_perm)**2)))

                    Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
                    soil.Lw = Lw_opt

                if ((soil.Lw != np.nan) & (soil.water_init != np.nan) & (soil.bulk_perm_init != np.nan)):
                    soil.df['bulk_perm'] = [pm.wunderlich_p(soil.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], soil.Lw) if (min(water_range) <= soil.water[x] <= max(water_range) ) else soil.bulk_perm[x] for x in range(soil.n_states)]
                    return soil.df.bulk_perm


            ### La freq no cambia
            elif (((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm <= 30e6)).all() or ((soil.df.frequency_perm > 100e6) & (soil.df.frequency_perm < 200e6)).all()):
                bulk_perm_inf(soil)                  # Celeste
                #TODO bulk_cond(soil)

                soil.df.loc[(np.isnan(soil.df['bulk_perm']) == True) & (((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm < 30e6)) or (soil.df.frequency_perm > 100e6) & (soil.df.frequency_perm < 200e6)), ['bulk_perm']] = pm.longmire_smith(soil.df.bulk_cond.values, soil.df.bulk_perm_inf.values, soil.df.frequency_perm.values)
                soil.df['bulk_perm'] = [pm.wunderlich_p(soil.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], soil.Lw) if (min(water_range) <= soil.water[x] <= max(water_range) ) else soil.bulk_perm[x] for x in range(soil.n_states)]
                return soil.df.bulk_perm.values


            elif (((soil.df.frequency_perm) >= 30e6) & (soil.df.frequency_perm <= 100e6)).any:  # Go over if any value is missing
                particle_density(soil)               # Yellow TODO revisar chart flow            
                air_perm(soil)                       # Celeste TODO revisar chart flow     
                solid_perm(soil)                     # Celeste TODO revisar chart flow 
                water_perm(soil)                     # Yellow
                CEC(soil)                            # Lila TODO ml approach

                soil.df.loc[(np.isnan(soil.df['bulk_perm']) == True) & (soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm <= 100e6), ['bulk_perm']] = pm.roth_mv(soil.water, soil.df.bulk_density.values, soil.df.particle_density.values, soil.df.air_perm.values, soil.df.solid_perm.values, soil.df.water_perm.values, soil.df.CEC.values)
                return soil.df.bulk_perm.values


            elif (((soil.df.frequency_perm) >= 200e6) & (soil.df.frequency_perm <= 14e9)).any:  # Go over if any value is missing
                particle_density(soil)               # Yellow TODO revisar chart flow            
                air_perm(soil)                       # Celeste TODO revisar chart flow     
                solid_perm(soil)                     # Celeste TODO revisar chart flow 
                water_perm(soil)                     # Yellow
                texture(soil)                        # Lila

                soil.df.loc[(np.isnan(soil.df['bulk_perm']) == True) & (soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 14e9), ['bulk_perm']] = pm.roth_wn(soil.water, soil.df.bulk_density, soil.df.particle_density, soil.df.air_perm, soil.df.solid_perm, soil.df.water_perm, soil.df.clay) 
                return soil.df.bulk_perm.values