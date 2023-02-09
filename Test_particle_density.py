# Test particle_density predictions

import pedophysics
import numpy as np
from predict.frequency_perm import frequency_perm
from predict.particle_density import particle_density


sample2 = pedophysics.Soil(water = np.append(np.random.rand(15)*40, np.full((4,1),np.nan)), 
                            bulk_perm = np.array([15, np.nan, 20, np.nan, 40]), 
                            instrument='TDR', 
                            particle_density=np.array([2, 2.2, 3, np.nan, 2.6]))

print('sample2.particle_density', sample2.particle_density, type(sample2.particle_density))
print('particle_density(sample2)', particle_density(sample2), type(particle_density(sample2)))


sample3 = pedophysics.Soil(water = np.append(np.random.rand(12)*40, np.full((4,1),np.nan)), 
                            sand = np.array([ np.nan, 30, 30,     20,     np.nan,  30,     np.nan, 20     ]), 
                            silt = np.array([ 10    , 30, np.nan, np.nan, np.nan,  np.nan, 20,     20     ]), 
                            clay = np.array([ np.nan, 30, 30,     np.nan, 20,      30,     20,     np.nan ]), 
                            orgm = np.array([ np.nan, 1,  np.nan, 1,      np.nan,  1,      0.5,    np.nan ]), 
                            particle_density=[2,      2,  2.2,    np.nan, np.nan,  np.nan, np.nan, np.nan] )


print('sample3.particle_density', sample3.particle_density, type(sample3.particle_density))
print('sample3.df.sand', sample3.df.sand)
print('sample3.df.silt', sample3.df.silt)
print('sample3.df.clay', sample3.df.clay)
print('particle_density(sample3)', particle_density(sample3), type(particle_density(sample3)))
print('sample3.df.sand', sample3.df.sand)
print('sample3.df.silt', sample3.df.silt)
print('sample3.df.clay', sample3.df.clay)

sample3_freq = frequency_perm(sample3)
print('sample3_freq', sample3_freq)
 
#print('sample2.df', sample2.df, len(sample2.df))
print(sample2.particle_density==sample3.particle_density)
print(sample2==sample3)