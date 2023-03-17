import numpy as np

def CEC(soil):
    """ 
        Cation exchange capacity prediction
    """
    #if np.isnan(soil.CEC == True).any:

    #    if (soil.orgm != np.nan) & (soil.OC == np.nan):
    #        soil.set_OC(ptf.orgm2oc(soil.orgm))

    #    if ((soil.sand == np.nan) or (soil.silt == np.nan) or (soil.clay == np.nan)) and (soil.texture != np.nan):
    #        print('Translating texture to fractions')
    #        ptf.texture_to_fraction(soil)
        
    #    if (soil.clay != np.nan) & (soil.OC != np.nan) & (soil.silt != np.nan) & (soil.depth != np.nan):
    #        return ptf.world_wosis_cecph7(soil.silt, soil.clay, soil.OC, soil.depth)

    #    else:
    #        print('CEC cannot be calculated. Please provide soil clay, silt, OC and depth')
    
    
    
    #soil.df['CEC'] = [10 for x in range(soil.n_states)]
    return soil.df.CEC.values