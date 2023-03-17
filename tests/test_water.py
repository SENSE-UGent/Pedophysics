import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'F:\\Users\\gamendoz\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import Water


print("################## Example1 ####################")     
                         #                      0     1     2     3       4       5       6       7
sample1 = Soil(water = np.array([   0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'TDR')

water1 = Water(sample1)
print('water1', water1)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1')
ax.plot(sample1.bulk_perm, sample1.water, 'bo',  markersize=4)
ax.plot(sample1.bulk_perm, water1, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)
plt.show()

print("################## Example1b ####################") 
                         #                      0     1     2     3       4       5       6       7
sample1b = Soil(water = np.array([  0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                    frequency_perm=np.array([   50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6]))

water1b = Water(sample1b)
print('water1b', water1b)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

print("################## Example1c ####################")
                         #                      0     1     2     3       4       5       6       7
sample1c = Soil(water = np.array([  0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_cond = 0.05,
                    frequency_perm=np.array([   50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6]))

water1c = Water(sample1c)
print('water1c', water1c)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1c')
ax.plot(sample1c.bulk_perm, sample1c.water, 'bo',  markersize=4)
ax.plot(sample1c.bulk_perm, sample1c.df.water.values, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)
plt.show()

print("################## Example2 ####################") # now when permittivity freq is not passed then return
                         #                      0    1    2   3    4  5  6   7   8  9   10
sample2 = Soil(water = np.array([  .20,  .31, 0.36, 0.38, 0.05]), 
                            bulk_perm=np.array([10,   15, 20,   25,   7, 1, 12, 22, 5, 20, 30   ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5)

water2 = Water(sample2)
print("water2", water2)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example2')
ax.plot(sample2.bulk_perm, sample2.water, 'bo',  markersize=4)
ax.plot(sample2.bulk_perm, water2, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example3 ####################")   # It fills all the missing water states so wat from cond is not called
                         #                      0    1    2   3    4  5  6   7   8  9   10
sample3 = Soil(water = np.array([   .20,  0.30, 0.35                                 ]), 
                            bulk_perm=np.array([10,   15,   20,   8.5, 8, 1, 12, 22, 5, 20, 30   ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'GPR')

water3 = Water(sample3)
print("water3", water3)


print("################## Example4 ####################")     
                         #                      0    1    2   3    4  5  6   7   8  9   10
sample4 = Soil(water = np.array([   .20, .30, 0.35]), 
                            bulk_perm=np.array([10,   15, 20, 8.5, 8, 1, 12, 22, 5, 20, 30   ]), 
                            bulk_density=1.7,
                            clay = 40,
                            solid_perm = 5,
                            instrument = 'TDR')

water4 = Water(sample4)
print("water4", water4)


print("################## Example5 ####################")     
                         #                      0     1   2   3    4  5  6   7   8  9   10
sample5 = Soil(water = np.array([   .20, .30, 0.35]), 
                            bulk_perm=np.array([10,   15, 20, 8.5, 8, 1, 12, 22, 5, 20, 30   ]), 
                            bulk_density=1.7,
                            alpha = 0.1,
                            frequency_perm= np.array([60e6]))

water5 = Water(sample5)
print("water5", water5)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example5')
ax.plot(sample5.bulk_perm, sample5.water, 'bo',            markersize=4)
ax.plot(sample5.bulk_perm, water5,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example6 ####################") 
                         #                      0     1   2     3     4  5  6   7   8  9   10
sample6 = Soil( bulk_perm=np.array([10,   15, 20,   25,   7, 1, 12, 22, 5, 20, 30   ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'GPR')

water6 = Water(sample6)
print("water6", water6)

print("################## Example6b ####################")
                         #                      0       1      2      3      4      5     6      7      8     9      10
sample6b = Soil( bulk_perm=np.array([10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                      frequency_perm = np.array([1e6,   2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  50e6, 100e6, 200e6]))

water6b = Water(sample6b)
print("water6b", water6b)

print("################## Example6c ####################")
                         #                      0       1      2      3      4      5     6      7      8     9      10
sample6c = Soil( bulk_perm=np.array([10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_cond = 0.1,
                      frequency_perm = np.array([1e6,   2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  50e6, 100e6, 200e6]))

water6c = Water(sample6c)
print("water6c", water6c)


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example6c')
ax.plot(sample6c.bulk_perm, water6c, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()