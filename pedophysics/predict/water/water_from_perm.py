from ..frequency_perm import FrequencyPerm
from ..bulk_perm_inf import BulkPermInf
from ..particle_density import ParticleDensity
from ..air_perm import AirPerm
from ..solid_perm import SolidPerm
from ..water_perm import WaterPerm
from ..cation_exchange_capacity import CEC
from ..texture import Texture

from pedophysics.pedophysical_models.water import RothCRIM, RothW, RothMV
from pedophysics.pedophysical_models.bulk_perm import WunderlichP, LongmireSmithP

import numpy as np
import warnings


######################################  Changing frequency module #####################################
def changing_freq(soil):    
    print("Water from perm Changing frequency module")
    BulkPermInf(soil)
    
    cond = []
    iterative_cond_cfr = np.arange(1e-6,  1,  1/(100000))
    for x in range(soil.n_states):

        if ~np.isnan(soil.bulk_perm[x]):
            LS_perm = LongmireSmithP(iterative_cond_cfr, soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x])           
            LS_eval = list(abs(LS_perm - [soil.bulk_perm[x]]*len(iterative_cond_cfr))) 
            cond_opt = iterative_cond_cfr[LS_eval.index(min(LS_eval))]
            cond.append(round(cond_opt, soil.roundn+3))

        else:
            cond.append(np.nan)
    soil.df['bulk_cond'] = [cond[x] if np.isnan(soil.df.bulk_cond[x]) else soil.df.bulk_cond[x] for x in range(soil.n_states)]
    print('soil.df.bulk_cond', soil.df.bulk_cond)


######################################  Fitting routine (fixed frequency) #####################################
def fitting(soil):
    print("Water from perm Fitting routine (fixed frequency)")
    WaterPerm(soil)                        # Yellow

    if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_perm_init)):
        soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])  
        soil.bulk_perm_init = min([soil.df.bulk_perm[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
        bulk_perm_final = max([soil.df.bulk_perm[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)])
        soil.bulk_perm_range = [soil.bulk_perm_init - (bulk_perm_final-soil.bulk_perm_init)/5, bulk_perm_final + (bulk_perm_final-soil.bulk_perm_init)/5]
        print('soil.bulk_perm_range', soil.bulk_perm_range)
        Lw_ = np.arange(-0.2, 0.8, 1/100)
        Wund_RMSE = []

        for i in range(len(Lw_)):
            wund_eval = [WunderlichP(soil.df.water[x], soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], Lw_[i]) if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else np.nan for x in range(soil.n_states)]
            Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_perm)**2)))
        Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
        soil.Lw = Lw_opt

    if ~np.isnan(soil.Lw) & ~np.isnan(soil.water_init) & ~np.isnan(soil.bulk_perm_init):
        Wundw = []
        iterative_wat = np.arange(0,  .65, 0.65/(4*65))
        for x in range(soil.n_states):
            
            if (min(soil.bulk_perm_range) <= soil.df.bulk_perm[x] <= max(soil.bulk_perm_range)) & ~np.isnan(soil.df.bulk_perm[x]):
                wundw_eval = list(abs(WunderlichP(iterative_wat, soil.bulk_perm_init, soil.water_init, soil.df.water_perm[x], soil.Lw) - [soil.df.bulk_perm[x]]*len(iterative_wat)))
                wat_opt = iterative_wat[wundw_eval.index(min(wundw_eval))]
                Wundw.append(round(wat_opt, soil.roundn))

            else:
                Wundw.append(np.nan)
        soil.df['water'] = [Wundw[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
        print('soil.df.water', soil.df.water)


######################################  non-fitting routine (fixed frequency) #####################################
def non_fitting(soil):
    print("Water from perm non-fitting routine (fixed frequency)")
    ParticleDensity(soil)               # Yellow TODO revisar chart flow            
    AirPerm(soil)                       # Celeste TODO revisar chart flow     
    SolidPerm(soil)                     # Celeste TODO revisar chart flow 
    WaterPerm(soil)                     # Yellow
    Texture(soil)                        # Lila
    CEC(soil)                            # Lila TODO ml approach
    BulkPermInf(soil)

    if ((soil.df.frequency_perm >= 100) & (soil.df.frequency_perm < 30e6)).all():
        cond = []
        iterative_cond_nf = np.arange(1e-6,  1,  1/100000)

        for x in range(soil.n_states):

            if ~np.isnan(soil.df.bulk_perm[x]):
                LS_perm = LongmireSmithP(iterative_cond_nf, soil.df.bulk_perm_inf[x], soil.df.frequency_perm[x])           
                LS_eval = list(abs(LS_perm - [soil.bulk_perm[x]]*len(iterative_cond_nf))) 
                cond_opt = iterative_cond_nf[LS_eval.index(min(LS_eval))]
                cond.append(round(cond_opt, soil.roundn+3))

            else:
                cond.append(np.nan)

        soil.df['bulk_cond'] = [cond[x] if np.isnan(soil.df.bulk_cond[x]) else soil.df.bulk_cond[x] for x in range(soil.n_states)]

    elif ((soil.df.frequency_perm >= 30e6) & (soil.df.frequency_perm < 100e6)).all():
        soil.df['water'] = [round(RothMV(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.CEC[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    elif ((soil.df.frequency_perm >= 100e6) & (soil.df.frequency_perm < 200e6)).all():
        soil.df['water'] = [round(RothCRIM(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]

    elif ( ((soil.df.frequency_perm >= 200e6) & (soil.df.frequency_perm <= 14e9))).all():  # Go over if any value is missing
        soil.df['water'] = [round(RothW(soil.df.bulk_perm[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.air_perm, soil.df.solid_perm[x], soil.df.water_perm[x], soil.df.clay[x]), soil.roundn) if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)] 


################################ Fixed frequnecy ##########################################
def fixed_freq(soil):
    print('Water from perm Fixed frequnecy')
            
    if (sum([1 if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_perm[x])) ) else 0 for x in range(soil.n_states)]) >= 3):
        fitting(soil)

    if (sum([1 if ( (np.isnan(soil.df.water[x])) & ((soil.df.frequency_perm[x] >= 100) & (soil.df.frequency_perm[x] <= 14e9)) ) else 0 for x in range(soil.n_states)]) >= 1):
        non_fitting(soil)


############################# Water prediction from permittivity ###################################### 
def WaterFromPerm(soil):
    """ 
        Soil volumetric water content prediction from real bulk relative dielectric permittivity
    """
    print("Water prediction from permittivity")
    FrequencyPerm(soil)

    if (np.isnan(soil.df.frequency_perm)).all():
        warnings.warn("All stages of 'frequency_perm' are Nan. Soil.water cannot be predicted by Soil.bulk_perm. Please define Soil.frequency_perm or Soil.instrument")
        return 

    elif np.all(soil.df.frequency_perm == soil.df.frequency_perm[0]):
        fixed_freq(soil)

    else :
        changing_freq(soil)