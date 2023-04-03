import numpy as np

from .water_from_cond import WaterFromCond
from .water_from_perm import WaterFromPerm

def Water(soil):
    """ 
        Soil volumetric water content prediction
    """
    # These conditions do not consider 
    if (sum([1 if np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_perm[x]) else 0 for x in range(soil.n_states)]) >= 1) :
        WaterFromPerm(soil)

    if (sum([1 if np.isnan(soil.df.water[x]) & ~np.isnan(soil.df.bulk_cond[x]) else 0 for x in range(soil.n_states)]) >= 1) :
        WaterFromCond(soil)        

    return soil.df.water.values

