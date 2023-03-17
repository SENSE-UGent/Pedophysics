from .water_perm import *
from .frequency_perm import *
from .bulk_cond import *
from .bulk_perm_inf import *
from .particle_density import *
from .air_perm import *
from .solid_perm import *
from .cation_exchange_capacity import *
from pedophysics.pedophysical_models.bulk_perm import *

import numpy as np


def BulkPerm(soil):
    """ 
        Real bulk dielectric permittivity prediction
    """
    if (np.isnan(soil.df.bulk_perm)).any():  # Go over if any value is missing
        FrequencyPerm(soil)

        if (np.isnan(soil.df.frequency_perm)).all():
            return

        elif (np.all(soil.df.frequency_perm == soil.df.frequency_perm[0])): # No changes in frequency

            if (sum([1 if ( ~np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_perm[x]) ) else 0 for x in range(soil.n_states)]) >= 3) :
                ################################################### Fitting routine ##################################################
                WaterPerm(soil)                        # Yellow

                if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_perm_init)):

                    soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
                    print('soil.water_init ', soil.water_init )  
                    soil.bulk_perm_init = min([soil.df.bulk_perm[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
                    water_final = max([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
                    soil.water_range = [soil.water_init - (water_final-soil.water_init)/5, water_final + (water_final-soil.water_init)/5]

                    Lw_ = np.arange(-0.2,  1,  (1.2) / 100 )
                    Wund_RMSE = []

                    for i in range(len(Lw_)):
                        wund_eval = [WunderlichP(soil.df.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], Lw_[i]) if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)]
                        Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm)**2)))

                    Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
                    soil.Lw = Lw_opt
                    print('Lw_opt', Lw_opt)

                if ((soil.Lw != np.nan) & (soil.water_init != np.nan) & (soil.bulk_perm_init != np.nan) & (soil.water_range != np.nan)):

                    soil.df['bulk_perm'] = [round(WunderlichP(soil.df.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], soil.Lw), soil.roundn) if (min(soil.water_range) <= soil.water[x] <= max(soil.water_range) ) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

            if ( (np.isnan(soil.df.bulk_perm)).any() & (((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm < 30e6)).all())):
            ########################## Non-fitting fixed frequency routine ############################
                print('entro low perm')
                BulkPermInf(soil)                  # Celeste
                BulkCond(soil)
                # Bulk density!!! TODO incluir porosity
                soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_cond[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) if (np.isnan(soil.df.bulk_perm[x])) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

            elif (np.isnan(soil.df.bulk_perm)).any() & ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm <= 14e9)).all(): 
                ParticleDensity(soil)               # Yellow TODO revisar chart flow            
                AirPerm(soil)                       # Celeste TODO revisar chart flow     
                SolidPerm(soil)                     # Celeste TODO revisar chart flow 
                WaterPerm(soil)                     # Yellow
                Texture(soil)                        # Lila
                CEC(soil)                            # Lila TODO ml approach

                if ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():

                    soil.df['bulk_perm'] = [round(RothMV(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]
                    # TODO integrating HydraProbe model for non CEC known states?

                elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():

                    soil.df['bulk_perm'] = [round(RothCRIM(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

                elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 14e9))).all():  # Go over if any value is missing

                    soil.df['bulk_perm'] = [round(RothW(soil.df.water[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.bulk_perm[x]) else soil.df.bulk_perm[x] for x in range(soil.n_states)] 

        else:
        ########################################### Changing frequency permittivity ############################################
            print('changing frequency permittivity')
            BulkPermInf(soil)                  # Celeste
            BulkCond(soil)
            
            soil.df['bulk_perm'] = [round(LongmireSmithP(soil.df.bulk_cond[x], soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]), soil.roundn) if (np.isnan(soil.df.bulk_perm[x])) else soil.df.bulk_perm[x] for x in range(soil.n_states)]

    return soil.df.bulk_perm.values
