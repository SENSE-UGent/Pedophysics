from pedophysics.predict import *

from pedophysics.pedophysical_models.bulk_cond import WunderlichC, LongmireSmithC, Fu

import numpy as np
import warnings


def non_fitting_dc(soil):
    print('Water from cond non-fitting routine, freq<=100Hz')
    wat = []
    wat_0 = 0
    wat_end = 0.65
    wat_step = 10**(-soil.roundn)
    iterative_wat= np.arange(wat_0,  wat_end, wat_step)

    for x in range(soil.n_states):
        Fu_diff = list(abs(Fu(iterative_wat, soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_cond[x], soil.df.solid_cond[x], soil.df.dry_cond[x], soil.df.sat_cond[x]) - np.array([soil.df.bulk_cond[x]]*len(iterative_wat))) if (soil.df.frequency_cond[x]<=100) else np.array([np.nan]*len(iterative_wat)) )

        if np.isnan(Fu_diff).all():
            wat.append(soil.df.water[x])
        
        else:
            wat_opt = iterative_wat[Fu_diff.index(np.nanmin(Fu_diff))]
            wat.append(wat_opt)

    soil.df['water'] = [round(wat[x], soil.roundn) if (np.isnan(soil.df.water[x]) & (soil.df.frequency_cond[x]<=100)) else soil.df.water[x] for x in range(soil.n_states) ]


def non_fitting_non_dc(soil):
    print('Water from cond non-fitting routine, freq>100Hz')
    cond_dc = []
    cond_dc_0 = 1e-6
    cond_dc_end = 1
    cond_dc_step = 10**(-(soil.roundn+3))
    iterative_cond_dc = np.arange(cond_dc_0, cond_dc_end, cond_dc_step)

    for x in range(soil.n_states):
        LS_eval = list(abs(LongmireSmithC(iterative_cond_dc, soil.df.frequency_cond[x]) - np.array([soil.df.bulk_cond[x]]*len(iterative_cond_dc))) if (soil.df.frequency_cond[x]>100) else np.array([np.nan]*len(iterative_cond_dc)))     

        if np.isnan(LS_eval).all():
            cond_dc.append(np.nan)

        else:
            cond_dc_opt = iterative_cond_dc[LS_eval.index(np.nanmin(LS_eval))]
            cond_dc.append(cond_dc_opt)

    print('cond_dc', cond_dc)
    ##########################################
    wat = []
    wat_0 = 0
    wat_end = 0.65
    wat_step = 10**(-soil.roundn)
    iterative_wat= np.arange(wat_0,  wat_end, wat_step)

    for x in range(soil.n_states):
        Fu_diff = list(abs(Fu(iterative_wat, soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_cond[x], soil.df.solid_cond[x], soil.df.dry_cond[x], soil.df.sat_cond[x]) - np.array([cond_dc[x]]*len(iterative_wat))) if (soil.df.frequency_cond[x]>100) else np.array([np.nan]*len(iterative_wat) ))
        
        if np.isnan(Fu_diff).all():
            wat.append(soil.df.water[x])
        
        else:
            wat_opt = iterative_wat[Fu_diff.index(np.nanmin(Fu_diff))]
            wat.append(wat_opt)

    soil.df['water'] = [round(wat[x], soil.roundn) if (np.isnan(soil.df.water[x]) & (soil.df.frequency_cond[x]>100)) else soil.df.water[x] for x in range(soil.n_states) ]


def non_fitting(soil):   
    print('Water from cond non-fitting routine') 
    Texture(soil)
    ParticleDensity(soil)
    WaterCond(soil)
    SolidCond(soil)

    if (sum([1 if np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) & (soil.df.frequency_cond[x] <= 100) else 0 for x in range(soil.n_states)]) >= 1):
        non_fitting_dc(soil)

    if (sum([1 if np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) & (soil.df.frequency_cond[x] > 100) else 0 for x in range(soil.n_states)]) >= 1):
        non_fitting_non_dc(soil)


######################################  Fitting routine (frequency_cond <= 100) #####################################
def fitting(soil):
    print('Water from cond Fitting routine (frequency_cond <= 100)')
    WaterCond(soil)                        # Yellow

    if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_cond_init)):
        soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
        soil.bulk_cond_init = min([soil.df.bulk_cond[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
        bulk_cond_final = max([soil.df.bulk_cond[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
        soil.bulk_cond_range = [soil.bulk_cond_init - (bulk_cond_final-soil.bulk_cond_init)/5, bulk_cond_final + (bulk_cond_final-soil.bulk_cond_init)/5]
        Lw_ = np.arange(-0.2, 0.8, 1/100)
        cond = []
        
        for i in range(len(Lw_)):
            wund_diff = np.absolute([WunderlichC(soil.df.water[x], soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], Lw_[i]) - soil.df.bulk_cond[x] if (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) and soil.df.frequency_cond[x]<100 else np.nan for x in range(soil.n_states)])    

            if np.isnan(wund_diff).all():
                cond.append(np.nan)

            else:
                cond.append(np.nansum(wund_diff))

        if np.isnan(cond).all():
            warnings.warn("Fitting cannot be done because because the pedophysical model 'WunderlichC' is badly evaluate. Please define Soil.water_cond", )
            return 

        else:
            Lw_opt = Lw_[cond.index(np.nanmin(cond))]
            soil.Lw = Lw_opt
            print('Lw_opt', Lw_opt)

    if ~np.isnan(soil.Lw) & ~np.isnan(soil.water_init) & ~np.isnan(soil.bulk_cond_init):
        Wundw = []
        iterative_wat = np.arange(0,  .65, 0.65/(4*65))
        for x in range(soil.n_states):

            if (min(soil.bulk_cond_range) <= soil.df.bulk_cond[x] <= max(soil.bulk_cond_range)) & ~np.isnan(soil.df.bulk_cond[x]):
                wundw_eval = list(abs(WunderlichC(iterative_wat, soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], soil.Lw) - [soil.df.bulk_cond[x]]*len(iterative_wat)))
                wat_opt = iterative_wat[wundw_eval.index(min(wundw_eval))]
                Wundw.append(round(wat_opt, soil.roundn))

            else:
                Wundw.append(np.nan)

        soil.df['water'] = [Wundw[x] if np.isnan(soil.df.water[x]) else soil.df.water[x] for x in range(soil.n_states)]
        print('soil.df.water (after fitting)', soil.df.water)


def WaterFromCond(soil):
    """ 
        Soil volumetric water content prediction from real bulk electrical conductivity
    """
    print("water from cond")
    FrequencyCond(soil)

    if (sum([1 if  ~np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) and soil.df.frequency_cond[x]<100 else 0 for x in range(soil.n_states)]) >= 3) :
        fitting(soil)
    
    if (sum([1 if np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) else 0 for x in range(soil.n_states)]) >= 1): #TODO else. sum ==2 is not considered
        non_fitting(soil)