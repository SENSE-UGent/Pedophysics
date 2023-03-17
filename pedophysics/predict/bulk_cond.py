import numpy as np
from pedophysics.pedophysical_models.bulk_cond import *
# from pedophysics.predict import FrequencyCond, ParticleDensity, WaterCond, Texture, SolidCond
import pedophysics

def BulkCond(soil):
    """ 
        Real bulk electrical conductivity prediction
    """
    if (np.isnan(soil.df.bulk_cond)).any():  # Go over if any value is missing
        FrequencyCond(soil)

        if (sum([1 if  ~np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) and soil.df.frequency_cond[x]<100 else 0 for x in range(soil.n_states)]) >= 3) :
            ################################################### Fitting routine ##################################################
            print('Cond fitting routine')
            WaterCond(soil)                        # Yellow

            if (np.isnan(soil.Lw) or np.isnan(soil.water_init) or np.isnan(soil.bulk_cond_init)):
                soil.water_init = min([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
                soil.bulk_cond_init = min([soil.df.bulk_cond[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
                water_final = max([soil.df.water[x] if ( (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) ) else np.nan for x in range(soil.n_states)])
                soil.water_range = [soil.water_init - (water_final-soil.water_init)/5, water_final + (water_final-soil.water_init)/5]
                Lw_ = np.arange(-0.2,  0.8,  1 / 100 )
                Wund_RMSE = []
                
                for i in range(len(Lw_)):
                    wund_eval = [WunderlichC(soil.df.water[x], soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], Lw_[i]) if (~np.isnan(soil.df.water[x])) & (~np.isnan(soil.df.bulk_cond[x])) and soil.df.frequency_cond[x]<100 else np.nan for x in range(soil.n_states)]
                    Wund_RMSE.append(np.sqrt(np.nanmean((np.array(wund_eval) - soil.df.bulk_cond)**2)))

                Lw_opt = Lw_[Wund_RMSE.index(min(Wund_RMSE))]
                soil.Lw = Lw_opt

            if ((soil.Lw != np.nan) & (soil.water_init != np.nan) & (soil.bulk_cond_init != np.nan) & (soil.water_range != np.nan)):
                soil.df['bulk_cond'] = [round(WunderlichC(soil.df.water[x], soil.bulk_cond_init, soil.water_init, soil.df.water_cond[x], soil.Lw), soil.roundn+3) if (min(soil.water_range) <= soil.water[x] <= max(soil.water_range)) and soil.df.frequency_cond[x]<100 else soil.df.bulk_cond[x] for x in range(soil.n_states)]

        if (np.isnan(soil.df.bulk_cond)).any():
        ########################################### Changing frequency conductivity ############################################
            print('no-fitting')
            Texture(soil)
            ParticleDensity(soil)
            WaterCond(soil)
            SolidCond(soil)

            bulk_cond_dc = [Fu(soil.df.water[x], soil.df.clay[x], soil.df.bulk_density[x], soil.df.particle_density[x], soil.df.water_cond[x], soil.df.solid_cond[x], soil.df.dry_cond[x], soil.df.sat_cond[x]) if (np.isnan(soil.df.bulk_cond[x])) else soil.df.bulk_cond[x] for x in range(soil.n_states)]
            soil.df['bulk_cond'] = [round(bulk_cond_dc[x], soil.roundn+3) if soil.df.frequency_cond[x]<100 else soil.df.bulk_cond[x] for x in range(soil.n_states)]
            soil.df['bulk_cond'] = [round(LongmireSmithC(bulk_cond_dc[x], soil.df.frequency_cond[x]), soil.roundn+3) if soil.df.frequency_cond[x]>=100 else soil.df.bulk_cond[x] for x in range(soil.n_states)]

    return soil.df.bulk_cond.values
    