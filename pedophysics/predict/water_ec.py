import numpy as np
from scipy.optimize import minimize, Bounds

from pedophysics.pedophysical_models.water_ec import SenGoode
from pedophysics.pedophysical_models.bulk_ec import Fu, Rhoades
from pedophysics.pedophysical_models.bulk_perm import Hilhorst
from pedophysics.utils.stats import R2_score

from .temperature import Temperature
from .particle_density import ParticleDensity
from .solid_ec import SolidEC
from .texture import Texture
from .frequency_ec import FrequencyEC
from .water_perm import WaterPerm


def WaterEC(soil):
    '''
    '''
    Temperature(soil)
    FrequencyEC(soil)

    # Condition for non-fitting approach
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.salinity[x]) for x in range(soil.n_states)) or sum(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)) == 1:
        non_fitting(soil)

    # Condition for fitting approach
    if sum(not np.isnan(soil.bulk_ec[x]) and (not np.isnan(soil.water[x]) or not np.isnan(soil.bulk_perm[x])) and np.isnan(soil.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting(soil)

    return soil.df.water_ec.values


######################################    non-fitting    ################################

def non_fitting(soil):
    """
    """
    # Condition for non-fitting approach using salinity
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.salinity[x]) for x in range(soil.n_states)):
        from_salinity(soil)

    # Condition for non-fitting approach using bulk_ec
    if any(np.isnan(soil.df.water_ec[x]) and not np.isnan(soil.df.water[x]) and not np.isnan(soil.df.bulk_ec[x]) for x in range(soil.n_states)):
        from_bulk_ec(soil)


######################################    from Salinity   #####################################

def from_salinity(soil):
    """
    """
    # Calculating and saving water_ec and its info
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated using SenGood function in predict.water_ec.from_salinity" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated using SenGood function in predict.water_ec.from_salinity"
                                else soil.info.water_ec[x] for x in range(soil.n_states)]
 
    soil.df['water_ec'] = [SenGoode(soil.df.temperature[x], soil.df.salinity[x]) if np.isnan(soil.df.water_ec[x]) 
                                else soil.df.water_ec[x] for x in range(soil.n_states)]


######################################    from Bulk EC   ######################################

def from_bulk_ec(soil):
    """
    """
    Texture(soil)
    ParticleDensity(soil)
    SolidEC(soil)

    # Defining minimization function to obtain water_ec
    def objective_wat_ec(water_ec, wat, clay, bulk_density, particle_density, solid_ec, dry_ec, sat_ec, EC):
        return abs(Fu(wat, clay, bulk_density, particle_density, water_ec, solid_ec, dry_ec, sat_ec) - EC)
    
    # Calculating optimal water_ec
    wat_ec = []
    for i in range(soil.n_states):
        res = minimize(objective_wat_ec, 0.14, args=(soil.df.water[i], soil.df.clay[i], soil.df.bulk_density[i], soil.df.particle_density[i], soil.df.solid_ec[i], 
                                                     soil.df.dry_ec[i], soil.df.sat_ec[i], soil.df.bulk_ec[i]), bounds=[(0, 2)] )
        wat_ec.append(np.nan if np.isnan(res.fun) else round(res.x[0], soil.roundn) )

    # Saving calculated water_ec and its info
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_ec.from_bulk_ec" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated using Fu function (reported R2=0.98) in predict.water_ec.from_bulk_ec"
                                 else soil.info.water_ec[x] for x in range(soil.n_states)]

    soil.df['water_ec'] = [round(wat_ec[x], soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water[x] for x in range(soil.n_states) ]


######################################   fitting from Bulk EC and (Bulk Perm or water)  #####################################

def fitting(soil):
    """
    """
    # Condition for fitting approach using Rhoades function
    if sum(not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting_rhoades(soil)
    
    # Condition for fitting approach using Rhoades function
    elif sum(not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 and np.isnan(soil.df.water_ec[x]) for x in range(soil.n_states)) >= 2:
        fitting_hilhorst(soil)


def fitting_rhoades(soil):
    """
    """
    # Selecting calibration data
    arg_EC_wn = np.array([soil.df.bulk_ec[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    arg_water_wn = np.array([soil.df.water[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.water[x]) else np.nan for x in range(soil.n_states)])
    
    # Removing NaNs from calibration data
    valid_indices = ~np.isnan(arg_EC_wn) & ~np.isnan(arg_water_wn)
    arg_EC = arg_EC_wn[valid_indices]
    arg_water = arg_water_wn[valid_indices]
    
    # Define the initial guesses
    bounds = Bounds([0.00001, 0], [2, 0.1])
    initial_guess_watec = 0.15
    initial_guess_s_ec = 0
    initial_guess_E = 1
    initial_guess_F = 0.38

    # Defining minimization function to obtain water_ec and s_ec while fixing E and F
    def objective_water_ec(params, wat, bulk_ec, E, F):
        water_ec, s_ec = params
        residuals = (Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec)**2
        return np.sum(residuals)

    # Calculating optimal water_ec and s_ec
    res1 = minimize(objective_water_ec, [initial_guess_watec, initial_guess_s_ec], args=(arg_water, arg_EC, initial_guess_E, initial_guess_F), bounds=bounds)
    best_water_ec, best_s_ecs = res1.x

    # Saving calculated s_ec and its info
    soil.info['s_ec'] = [str(soil.info.s_ec[x]) + "--> Calculated by fitting Rhoades function in predict.water_ec.fitting_rhoades" if np.isnan(soil.df.s_ec[x])
                            or soil.info.s_ec[x] == str(soil.info.s_ec[x]) + "--> Calculated by fitting Rhoades function in predict.water_ec.fitting_rhoades"
                            else soil.info.s_ec[x] for x in range(soil.n_states)]
    
    soil.df['s_ec'] = [round(best_s_ecs, soil.roundn+3) if np.isnan(soil.df.s_ec[x]) else soil.df.s_ec[x] for x in range(soil.n_states) ]

    # Defining minimization function to obtain E and F while fixing water_ec and s_ec
    def objective_others(params, wat, bulk_ec, water_ec, s_ec):
        E, F = params
        residuals = np.sum((Rhoades(wat, water_ec, s_ec, E, F) - bulk_ec)**2)
        return residuals

    # Calculating optimal E and F
    res2 = minimize(objective_others, [initial_guess_E, initial_guess_F], args=(arg_water, arg_EC, best_water_ec, best_s_ecs))
    best_E, best_F = res2.x
    soil.E = best_E
    soil.F = best_F

    # Calculating the R2 score of the fitting
    R2 = round(R2_score(arg_EC, Rhoades(arg_water, best_water_ec, best_s_ecs, best_E, best_F)), soil.roundn)
    
    # Saving calculated water_ec and its info with R2
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2 = "+str(R2)+") Rhoades function in predict.water_ec.fitting_rhoades" if np.isnan(soil.df.water_ec[x]) 
                            or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2 = "+str(R2)+") Rhoades function in predict.water_ec.fitting_rhoades"
                            else soil.info.water_ec[x] for x in range(soil.n_states)]

    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water_ec[x] for x in range(soil.n_states) ]


def fitting_hilhorst(soil):
    """
    """
    WaterPerm(soil)

    # Selecting calibration data
    arg_EC_wn = np.array([soil.df.bulk_ec[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                            else np.nan for x in range(soil.n_states)])
    arg_bulk_perm_wn = np.array([soil.df.bulk_perm[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                              else np.nan for x in range(soil.n_states)])
    arg_water_perm_wn = np.array([soil.df.water_perm[x] if not np.isnan(soil.df.bulk_ec[x]) and not np.isnan(soil.df.bulk_perm[x]) and soil.df.bulk_perm[x]>=10 
                               else np.nan for x in range(soil.n_states)])

    # Removing NaNs from calibration data
    valid_indices = ~np.isnan(arg_EC_wn) & ~np.isnan(arg_bulk_perm_wn)
    arg_EC = arg_EC_wn[valid_indices]
    arg_bulk_perm = arg_bulk_perm_wn[valid_indices]
    arg_water_perm = arg_water_perm_wn[valid_indices]
    
    # Define the initial guesses
    bounds = Bounds([0.00001, -10], [2, 10])
    initial_guess_offset_perm = 4
    initial_guess_watec = 0.15

    # Defining minimization function
    def objective_water_ec(param, bulk_perm, bulk_ec, water_perm):
        water_ec, offset_perm = param
        residuals = (Hilhorst(bulk_ec, water_ec, water_perm, offset_perm) - bulk_perm)**2
        return np.sum(residuals)

    # Calculating optimal water_ec and offset_perm
    res = minimize(objective_water_ec, [initial_guess_watec, initial_guess_offset_perm], args=(arg_bulk_perm, arg_EC, arg_water_perm), bounds=bounds)
    best_water_ec, best_offset_perm = res.x

    # Saving calculated offset_perm and its info
    soil.info['offset_perm'] = [str(soil.info.offset_perm[x]) + "--> Calculated by fitting Hilhorst function in predict.water_ec.fitting_hilhorst" if np.isnan(soil.df.offset_perm[x]) 
                                or soil.info.offset_perm[x] == str(soil.info.offset_perm[x]) + "--> Calculated by fitting Hilhorst function in predict.water_ec.fitting_hilhorst"
                                 else soil.info.offset_perm[x] for x in range(soil.n_states)]
    
    soil.df['offset_perm'] = [round(best_offset_perm, soil.roundn+3) if np.isnan(soil.df.offset_perm[x]) else soil.df.offset_perm[x] for x in range(soil.n_states) ]

    # Calculating the R2 score of the fitting
    R2 = round(R2_score(arg_bulk_perm, Hilhorst(arg_EC, best_water_ec, arg_water_perm, best_offset_perm)), soil.roundn)
    
    # Saving calculated water_ec and its info with R2
    soil.info['water_ec'] = [str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") Hilhorst function in predict.water_ec.fitting_hilhorst" if np.isnan(soil.df.water_ec[x]) 
                             or soil.info.water_ec[x] == str(soil.info.water_ec[x]) + "--> Calculated by fitting (R2="+str(R2)+") Hilhorst function in predict.water_ec.fitting_hilhorst"
                                 else soil.info.water_ec[x] for x in range(soil.n_states)]
    
    soil.df['water_ec'] = [round(best_water_ec, soil.roundn+3) if np.isnan(soil.df.water_ec[x]) else soil.df.water_ec[x] for x in range(soil.n_states) ]