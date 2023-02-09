# Water_permittivity test


import pedophysics
import numpy as np
from predict.frequency_perm import frequency_perm
from predict.water_perm import water_perm
import pandas as pd

pd.set_option("display.max_rows", None, "display.max_columns", None)

sample1 = pedophysics.Soil(water = np.array([0.2, 0.3, 0.5]), 
    water_perm = np.array([15,  np.nan, 20,  np.nan, 40,     15,  np.nan, 20,  np.nan, 40 ]), 
    instrument = 'HydraProbe', 
   temperature = np.array([220, np.nan, 215, np.nan, np.nan, 215, np.nan, 220, 230,    250, 250]),
frequency_perm = np.array([1e5, 1e5,    1e6, 1e6,    1e7,    1e7, 1e8,    1e8, 1e9,    1e9, np.nan, np.nan, np.nan]) )


print('predict.frequency_perm(sample1)', frequency_perm(sample1))                          
print('predict.water_perm(sample1)', water_perm(sample1))                          



