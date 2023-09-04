# Test instrument to frequency prediction

import numpy as np
import pandas as pd
from pedophysics.simulate_old import Soil
from pedophysics.predict import FrequencyPerm

# Checking frequency_perm module usage
sample1 = Soil(water = np.array([0.2, 0.3, 0.5]), 
                            bulk_perm = np.array([15, np.nan, 20, np.nan, 40]), 
                            instrument='TDR', 
                            frequency_perm=np.append(np.random.rand(10)*40, np.full((4,1),np.nan)))

print('sample1.frequency_perm', sample1.frequency_perm, len(sample1.frequency_perm))
print('sample1.instrument', sample1.instrument, len(sample1.instrument))
#print('sample1.df', sample1.df, len(sample1.df))
sample1_freq = FrequencyPerm(sample1)
print('sample1_freq', sample1_freq)
