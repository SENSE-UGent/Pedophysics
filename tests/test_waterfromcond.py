import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, 'F:\\Users\\gamendoz\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import Water

a = np.absolute(np.nan)
print('a', a)

print("################## Example1 ####################")     
                         #                      0     1     2     3       4       5       6       7
sample1 = Soil(water = np.array([               0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07]), 
                            bulk_cond=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5 ])*1e-3, 
                            water_cond = 0.05)

water1 = Water(sample1)
print('water1', water1)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1')
ax.plot(sample1.bulk_cond, sample1.water, 'bo',  markersize=4)
ax.plot(sample1.bulk_cond, water1, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1*1e-3, 32*1e-3)
ax.set_ylim(-0.1, 0.65)
plt.show()

print("################## Example1b ####################") 
                         #                      0     1     2     3       4       5       6       7     8       9
sample1b = Soil(water = np.array([              0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_cond=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5
                            )

water1b = Water(sample1b)
print('water1b', water1b)


print("################## Example1c ####################")
                        #                       0     1     2     3       4       5       6       7     8       9
sample1c = Soil(water = np.array([              0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_cond=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_cond = 0.05,
                    frequency_cond=np.array([   50e1, 50e2, 50e2, 200e2,  200e2,  200e2,  50e2,   50e1, 50e1,   200e2]))

water1c = Water(sample1c)
print('water1c', water1c)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1c')
ax.plot(sample1c.df.bulk_cond, sample1c.water, 'bo',  markersize=4)
ax.plot(sample1c.df.bulk_cond, sample1c.df.water.values, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()

print("################## Example2 ####################") 
                         #                      0     1     2     3     4     5    6     7     8    9   10
sample2 = Soil(water = np.array([               0.20, 0.31, 0.36, 0.38, 0.05                                 ]), 
                            bulk_cond=np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20, 30   ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand')

water2 = Water(sample2)
print("water2", water2)



print("################## Example3 ####################")  
                         #                      0     1     2     3    4  5  6   7   8  9   10
sample3 = Soil(water = np.array([              0.20, 0.30, 0.35                                 ]), 
                            bulk_cond=np.array([10,   15,   20,   8.5, 8, 1, 12, 22, 5, 20, 30   ])*1e-3, 
                            bulk_density=1.7,
                            water_cond=0.5,
                            texture = 'Sand')

water3 = Water(sample3)
print("water3", water3)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example3')
ax.plot(sample3.bulk_cond, sample3.water, 'bo',  markersize=4)
ax.plot(sample3.bulk_cond, water3, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)

plt.show()


print("################## Example4 ####################")     
                         #                      0      1     2   3    4  5  6   7   8  9   10
sample4 = Soil(water = np.array([               0.20,  0.30, 0.35]), 
                            bulk_cond=np.array([10,    15,   20, 8.5, 8, 1, 12, 22, 5, 20, 30   ])*1e-3, 
                            bulk_density=1.7,
                            clay = 40)

water4 = Water(sample4)
print("water4", water4)


print("################## Example5 ####################")     
                         #                      0     1     2   3    4  5  6   7   8  9   10
sample5 = Soil(water = np.array([               0.20, 0.30, 0.35]), 
                            bulk_cond=np.array([10,   15,   20, 8.5, 8, 1, 12, 22, 5, 20, 30   ])*1e-3, 
                            bulk_density=1.7,
                            clay = 40,
                            water_cond=0.4,
                            frequency_cond=5e3
                            )

water5 = Water(sample5)
print("water5", water5)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example5')
ax.plot(sample5.bulk_cond, sample5.water, 'bo',            markersize=4)
ax.plot(sample5.bulk_cond, water5,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)

plt.show()


print("################## Example6 ####################") 
                         #          0     1     2     3       4     5      6      7      8     9      10
sample6 = Soil( bulk_cond=np.array([10,   15,   20,   25,     7,    1,     12,    22,    5,    20,    30   ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
            water_cond = np.array([ 0.05, 0.06, 0.07, np.nan, 0.01, 0.1]))

water6 = Water(sample6)
print("water6", water6)

print("################## Example6b ####################")
                         #                       0     1      2      3      4      5     6      7      8     9      10
sample6b = Soil( bulk_cond=np.array([            10,   15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                              bulk_density=1.7,
                              texture = 'Sand',
                              water_cond = 0.1,
                        frequency_cond=np.array([1,    2,     2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

water6b = Water(sample6b)
print("water6b", water6b)

print("################## Example6c ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample6c = Soil( bulk_cond=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Clay',
                            water_cond = 0.1,
                    frequency_cond = np.array([ 1,    2,      2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

water6c = Water(sample6c)
print("water6c", water6c)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example6c')
ax.plot(sample6c.bulk_cond, water6c, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)

plt.show()


print("################## Example7 ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample7 = Soil( bulk_cond=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_cond = 0.1)

water7 = Water(sample7)
print("water7", water7)

print("################## Example7b ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample7b = Soil( bulk_cond=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_cond = 0.01,
                      frequency_cond = np.array([500]))

water7b = Water(sample7b)
print("water7b", water7b)

print("################## Example7c ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample7c = Soil( bulk_cond=np.array([           10,    0,    np.nan, np.nan,   7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_cond = 0.1,
                      frequency_cond = np.array([2e6]))

water7c = Water(sample7c)
print("water7c", water7c)

print("################## ExampleSandValthe ####################")
                         #                      0    1    2     3     4      5      6      7    8    9      
samplev = Soil( bulk_cond=np.array([            3,   8,   15,   20,   22,    7,     12,    18,  10,  2     ])*1e-3, 
                            bulk_density=1.4,
                            texture = 'Sand',
                            CEC = 1.6,
                            water_cond = 0.1,
                      frequency_cond = np.array([50, 500, 5000, 2000, 50000, 55000, 25000, 100, 10, 50]))

waterv = Water(samplev)
print("waterv", waterv)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('ExampleSandValthe')
ax.plot(samplev.bulk_cond, waterv, 'ro', alpha=0.3, markersize=8)

ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)

plt.show()