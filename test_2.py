# This test is for testing and finding an optimal code writing when aiming to evaluate functions

import pedophysics
import numpy as np
import pandas as pd
import pedotransfer_functions as ptf
import time

pd.set_option("display.max_rows", None, "display.max_columns", None)

def complete(silt, clay):
    return 100 - silt-clay

sample1 = pedophysics.Soil(water = np.array([0.2, 0.3, 0.31]), 
                            orgm = 1., 
                            instrument='TDR', 
                            silt=np.random.rand(10000)*40, 
                            clay = np.random.rand(6000)*50)
#print(sample1.sand)
#print('es?', [x for x in range(sample1.n_states) if (sample1.sand[x] != np.nan) ])

start = time.time()
sample1.df.sand = [lambda x: 100 - sample1.clay[x] - sample1.silt[x] if ((np.isnan(sample1.sand[x]) == True) & (np.isnan(sample1.silt[x]) == False) & (np.isnan(sample1.clay[x]) == False)) else np.nan for x in range(sample1.n_states)]
end = time.time()
print('1', end - start)
test1 = sample1.df.sand
#print('test1', test1)
sample1.df.sand = [0 for x in range(sample1.n_states)]
#print('test1', test1)

start = time.time()
t1 = [100 - sample1.clay[x] - sample1.silt[x] if ((np.isnan(sample1.sand[x]) == True) & (np.isnan(sample1.silt[x]) == False) & (np.isnan(sample1.clay[x]) == False)) else np.nan for x in range(sample1.n_states)]
sample1.df.sand = t1
end = time.time()
print('2', end - start)
test1 = sample1.df.sand
#print('test1', test1)
sample1.df.sand = [0 for x in range(sample1.n_states)]

start = time.time()
sample1.df.sand = [ptf.complete(sample1.clay[x], sample1.silt[x]) if ((np.isnan(sample1.sand[x]) == True) & (np.isnan(sample1.silt[x]) == False) & (np.isnan(sample1.clay[x]) == False)) else np.nan for x in range(sample1.n_states)]
end = time.time()
print('3', end - start)
test2 = sample1.df.sand
#print('test2', test2)
sample1.df.sand = [0 for x in range(sample1.n_states)]

start = time.time()
sample1.df['sand'] = [ptf.complete(sample1.clay[x], sample1.silt[x]) if ((np.isnan(sample1.sand[x]) == True) & (np.isnan(sample1.silt[x]) == False) & (np.isnan(sample1.clay[x]) == False)) else np.nan for x in range(sample1.n_states)]
end = time.time()
print('4', end - start)
test3 = sample1.df.sand
sample1.df.sand = [0 for x in range(sample1.n_states)]

start = time.time()
sample1.df['sand'] = [complete(sample1.clay[x], sample1.silt[x]) if ((np.isnan(sample1.sand[x]) == True) & (np.isnan(sample1.silt[x]) == False) & (np.isnan(sample1.clay[x]) == False)) else np.nan for x in range(sample1.n_states)]
end = time.time()
print('5', end - start)
test4 = sample1.df.sand
#print(sample1.df.sand) # hacete bien el boludo
sample1.df.sand = [0 for x in range(sample1.n_states)]

start = time.time()
sample1.df.loc[((np.isnan(sample1.df['sand']) == True) | (np.isnan(sample1.df['silt']) == False) |(np.isnan(sample1.df['clay']) == False)), ['sand']] = complete(sample1.clay, sample1.silt) 
end = time.time()
print('6', end - start)
sample1.df.sand = [0 for x in range(sample1.n_states)]

start = time.time()
sample1.df.loc[((np.isnan(sample1.df['sand']) == True) | (np.isnan(sample1.df['silt']) == False) |(np.isnan(sample1.df['clay']) == False)), ['sand']] = list(map(complete, sample1.clay, sample1.silt)) 
end = time.time()
print('7', end - start)
sample1.df.sand = [0 for x in range(sample1.n_states)]

print(test1.all() == test2.all() == test3.all() == test4.all())


#################################################################################



sample2 = pedophysics.Soil(water = np.append(np.random.rand(20)*40, np.full((4,1),np.nan)), 
                            sand = np.array([ np.nan, np.nan, 30,  20,     np.nan,  np.nan,  20,     20, 50,     10     ]), 
                            silt = np.array([ np.nan, np.nan, 30,  np.nan, 20,      30,      20,     20, np.nan, np.nan     ]), 
                            clay = np.array([ np.nan, 30,     30,  np.nan, np.nan,  30,      np.nan, 20, np.nan, np.nan     ]), 
                            orgm = np.array([ np.nan, 1,  np.nan, 1,      np.nan,  1,   0.5,    np.nan ]), 
                            particle_density=[2, 2,  2.2,    3])

start = time.time()
print('complete(sample2.clay, sample2.silt)', complete(sample2.clay, sample2.silt), len(complete(sample2.clay, sample2.silt)) )
print('sample2.df[sand]', sample2.df['sand'])
print('condition', ((np.isnan(sample2.df['sand']) == True) | (np.isnan(sample2.df['silt']) == False) |(np.isnan(sample2.df['clay']) == False)))
sample2.df.loc[((np.isnan(sample2.df['sand']) == True) | (np.isnan(sample2.df['silt']) == False) |(np.isnan(sample2.df['clay']) == False)), ['sand']] = complete(sample2.clay, sample2.silt)

end = time.time()
print('8', end - start)
sample2.df.sand
