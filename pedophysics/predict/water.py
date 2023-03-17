import numpy as np

from .water_perm import *
from .water_cond import *
from .frequency_perm import *
from .frequency_cond import *
from .bulk_cond import *
from .bulk_perm_inf import *
from .particle_density import *
from .air_perm import *
from .solid_perm import *
from .solid_cond import *
from .cation_exchange_capacity import *
from .texture import *

from pedophysics.pedophysical_models.water import *
from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP
from pedophysics.pedophysical_models.bulk_cond import WunderlichC, LongmireSmithC, Fu


def Water(soil):
    """ 
        Soil volumetric water content prediction
    """
    if (np.isnan(soil.df.water)).any():
        WaterFromPerm(soil)

    if (np.isnan(soil.df.water)).any():
        WaterFromCond(soil)        

    return soil.df.water.values


def WaterFromCond(soil):
    """ 
        Soil volumetric water content prediction from real bulk electrical conductivity
    """
    print("water from cond")
    FrequencyCond(soil)

    if (sum([1 if  ~np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) and soil.df.frequency_cond[x]<100 else 0 for x in range(soil.n_states)]) >= 3) :
    ##################### Fitting routine ###########################
        print('Water from cond fitting routine')
        WaterCond(soil)                        # Yellow

        if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_cond_init)):
            soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
            soil.bulk_cond_init = min([soil.df.bulk_cond[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
            water_final = max([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
            soil.bulk_cond_range = [soil.water_init - (water_final-soil.water_init)/5, water_final + (water_final-soil.water_init)/5]
            Lw_ = np.arange(-0.2, 0.8, 1/100)
            Wund_RMSE = []
            
            for i in range(len(Lw_)):
                wund_eval = [WunderlichC(soil.df.water[x], soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], Lw_[i]) if (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) and soil.df.frequency_cond[x]<100 else np.nan for x in range(soil.n_states)]
                Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_cond)**2)))

            Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
            soil.Lw = Lw_opt

        if ((soil.Lw != np.nan) & (soil.water_init != np.nan) & (soil.bulk_cond_init != np.nan) & (soil.bulk_cond_range != np.nan)):
            Wundw = []
            iterative_wat = np.arange(0,  .65, 0.65/(4*65))
            for x in range(soil.n_states):
                
                if (min(soil.bulk_cond_range) <= soil.bulk_cond[x] <= max(soil.bulk_cond_range)) & ~np.isnan(soil.bulk_cond[x]):
                    wundw_eval = [abs(WunderlichC(iterative_wat, soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], soil.Lw) - [soil.df.bulk_cond[x]]*len(iterative_wat))]
                    wat_opt = iterative_wat[wundw_eval.index(min(wundw_eval))]
                    Wundw.append(round(wat_opt, soil.roundn))

                else:
                    Wundw.append(np.nan)

            soil.df['water'] = [Wundw[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
            print('soil.df.water', soil.df.water)

    if (np.isnan(soil.df.water)).any():
        print('Water from cond non-fitting routine')
        Texture(soil)
        ParticleDensity(soil)
        WaterCond(soil)
        SolidCond(soil)

        if (soil.df.frequency_cond < 100).any():
        ######### DC frequency ########
            print('Water from cond non-fitting routine, freq<100Hz')
            wat = []
            iterative_wat= np.arange(0,  .65, 0.65/(4*65))

            for x in range(soil.n_states):
                fu_eval = [abs(Fu(iterative_wat, soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_cond[x], soil.df.solid_cond[x], soil.df.dry_cond[x], soil.df.sat_cond[x]) - [soil.df.bulk_cond[x]]*len(iterative_wat)) if soil.df.frequency_cond[x]<100 else np.nan]
                wat_opt = iterative_wat[fu_eval.index(min(fu_eval))]
                wat.append(wat_opt)

            soil.df['water'] = [round(wat[x], soil.roundn) if (np.isnan(soil.df.water[x]) & soil.df.frequency_cond[x]<100) else soil.df.water[x] for x in range(soil.n_states) ]

        if (soil.df.frequency_cond >= 100).any():
        ######## non-DC frequency #######
            print('Water from cond non-fitting routine, freq>100Hz')
            cond = []
            iterative_cond = np.arange(1e-6,  1,  1/100000)

            for x in range(soil.n_states):
                LS_eval = [abs(LongmireSmithC(iterative_cond, soil.df.frequency_cond[x]) - [soil.bulk_cond[x]]*len(iterative_cond)) if soil.df.frequency_cond[x]>=100 else np.nan]        
                cond_opt = iterative_cond[LS_eval.index(min(LS_eval))]
                cond.append(cond_opt)


            soil.df['water'] = [round(Fu(cond[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_cond[x], soil.df.solid_cond[x], soil.df.dry_cond[x], soil.df.sat_cond[x]), soil.roundn) if (np.isnan(soil.df.water[x])) else soil.df.water[x] for x in range(soil.n_states)]
            print('soil.df.water', soil.df.water)



######################################################################################################################################################


def WaterFromPerm(soil):
    """ 
        Soil volumetric water content prediction from real bulk relative dielectric permittivity
    """
    print("water from perm")
    FrequencyPerm(soil)

    if (np.isnan(soil.df.frequency_perm)).all():
        return ### TODO solve Warnings

    elif np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]): # No changes in frequency
        print('Water from perm No changes in frequency')

        if ((sum([1 if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else 0 for x in range(soil.n_states)]) >= 3)) :
        ########################### Fitting routine ###########################
            print("Water from perm fitting routine")
            WaterPerm(soil)                        # Yellow

            if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_perm_init)):
                soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])  
                soil.bulk_perm_init = min([soil.df.bulk_perm[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
                bulk_perm_final = max([soil.df.bulk_perm[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
                soil.bulk_perm_range = [soil.bulk_perm_init - (bulk_perm_final-soil.bulk_perm_init)/5, bulk_perm_final + (bulk_perm_final-soil.bulk_perm_init)/5]
                Lw_ = np.arange(-0.2, 0.8, 1/100)
                Wund_RMSE = []

                for i in range(len(Lw_)):
                    wund_eval = [WunderlichP(soil.df.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], Lw_[i]) if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)]
                    Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm)**2)))

                Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
                soil.Lw = Lw_opt

            if ((soil.Lw != np.nan) & (soil.water_init != np.nan) & (soil.bulk_perm_init != np.nan)):
                Wundw = []
                iterative_wat = np.arange(0,  .65, 0.65/(4*65))
                for x in range(soil.n_states):
                    
                    if (min(soil.bulk_perm_range) <= soil.bulk_perm[x] <= max(soil.bulk_perm_range)) & ~np.isnan(soil.bulk_perm[x]):
                        wundw_eval = list(abs(WunderlichP(iterative_wat, soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], soil.Lw) - [soil.bulk_perm[x]]*len(iterative_wat)))
                        wat_opt = iterative_wat[wundw_eval.index(min(wundw_eval))]
                        Wundw.append(round(wat_opt, soil.roundn))

                    else:
                        Wundw.append(np.nan)
                soil.df['water'] = [Wundw[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
                print('soil.df.water', soil.df.water)

        if (np.isnan(soil.df.water)).any() & ((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm <= 14e9)).all(): 
        #################################### Non- fitting routine ####################################
            print("Water from perm NON-fitting routine")
            ParticleDensity(soil)               # Yellow TODO revisar chart flow            
            AirPerm(soil)                       # Celeste TODO revisar chart flow     
            SolidPerm(soil)                     # Celeste TODO revisar chart flow 
            WaterPerm(soil)                     # Yellow
            Texture(soil)                        # Lila
            CEC(soil)                            # Lila TODO ml approach
            FrequencyPerm(soil)
            BulkPermInf(soil)

            if ((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm < 30e6)).all():
                cond = []
                iterative_cond_nf = np.arange(1e-6,  1,  1/100000)

                for x in range(soil.n_states):
                    LS_eval = list(abs(LongmireSmithP(iterative_cond_nf, soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x]) - [soil.bulk_perm[x]]*len(iterative_cond_nf)))        
                    cond_opt = iterative_cond_nf[LS_eval.index(min(LS_eval))]
                    cond.append(round(cond_opt, soil.roundn+3))

                soil.df['bulk_cond'] = [cond[x] if np.isnan(soil.df.bulk_cond[x]) else soil.df.bulk_cond[x] for x in range(soil.n_states)]
                print('soil.df.bulk_cond', soil.df.bulk_cond)

            elif ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():
                soil.df['water'] = [round(RothMV(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
                # TODO integrating HydraProbe model for non CEC known states?
                print('soil.df.water', soil.df.water)

            elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():
                soil.df['water'] = [round(RothCRIM(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
                print('soil.df.water', soil.df.water)

            elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 14e9))).all():  # Go over if any value is missing
                soil.df['water'] = [round(RothW(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)] 
                print('soil.df.water', soil.df.water)

    else :
    ##################### Changing frequency ###################
        print("Water from perm changing frequency")
        FrequencyPerm(soil)
        BulkPermInf(soil)
        
        cond = []
        iterative_cond_cfr = np.arange(1e-6,  1,  1/100000)

        for x in range(soil.n_states):
            LS_perm = LongmireSmithP(iterative_cond_cfr, soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x])
            LS_eval = list(abs(LS_perm - [soil.bulk_perm[x]]*len(iterative_cond_cfr)))        
            cond_opt = iterative_cond_cfr[LS_eval.index(min(LS_eval))]
            cond.append(round(cond_opt, soil.roundn+3))

        soil.df['bulk_cond'] = [cond[x] if np.isnan(soil.df.bulk_cond[x]) else soil.df.bulk_cond[x] for x in range(soil.n_states)]
        print('soil.df.bulk_cond', soil.df.bulk_cond)