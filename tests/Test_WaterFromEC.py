import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
#sys.path.insert(0, 'C:\\Users\\gasto\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import Water
from pedophysics.pedophysical_models.bulk_ec import WunderlichEC

##############
########################## Testing DC frequency, fitting and non-fitting ###############################
##############

print("################## Example1 ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample1WEC = Soil( bulk_ec=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_ec = 0.1)

water1 = Water(sample1WEC) #  
print("water1", water1) #    [0.264 0.335 0.394 0.447 0.213 0.065 0.293 0.394 0.173 0.394 0.416]
#                            
sample1WEC.info.to_excel('sample1WEC_info.xlsx')
sample1WEC.df.to_excel('sample1WEC_df.xlsx')


print("################## Example1b ####################") 
                         #          0     1     2     3       4     5      6      7      8     9      10
sample1bWEC = Soil( bulk_ec=np.array(  [10,   15,   20,   25,     7,    1,     12,    22,    5,    20,    30   ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
            water_ec = np.array( [ 0.05, 0.06, 0.07, np.nan, 0.01, 0.1]))

water1b = Water(sample1bWEC) #         [0.35  0.413 0.459   nan 0.442 0.065   nan   nan   nan   nan   nan]
print("water1b", water1b) #            
#                                     
sample1bWEC.info.to_excel('sample1bWEC_info.xlsx')
sample1bWEC.df.to_excel('sample1bWEC_df.xlsx')


print("################## Example2 ####################") 
                         #                    0     1     2     3     4     5    6     7     8    9     10
sample2WEC = Soil(water = np.array([          0.20, 0.31, 0.36, 0.38, 0.05                                   ]), 
                            bulk_ec=np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20,   30   ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand')

water2 = Water(sample2WEC) #               
print("water2", water2) #                    
#                                            [0.2   0.31  0.36  0.38  0.05  0.    0.186 0.396 0.    0.357 0.539]
sample2WEC.info.to_excel('sample2WEC_info.xlsx')
sample2WEC.df.to_excel('sample2WEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example2')
ax.plot(sample2WEC.bulk_ec, sample2WEC.water, 'bo',  markersize=4)
ax.plot(sample2WEC.bulk_ec, water2, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()


print("################## Example3 ####################")  
                         #                    0     1     2     3     4     5    6     7     8    9   10
sample3WEC = Soil(water = np.array([          0.20, 0.31, 0.36, 0.38, 0.05                                 ]), 
                            bulk_ec=np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20, 30   ])*1e-3, 
                            bulk_density=1.7,
                            water_ec=0.5,
                            texture = 'Sand')

water3 = Water(sample3WEC) # [0.2   0.31  0.36  0.38  0.05  0.    0.214 0.373 0.    0.349 0.446]
print("water3", water3) # 

sample3WEC.info.to_excel('sample3WEC_info.xlsx')
sample3WEC.df.to_excel('sample3WEC_df.xlsx')


print("################## Example4 ####################")     
                         #                      0     1     2     3       4       5       6       7
sample4WEC = Soil(water =                      [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07], 
                            bulk_ec=np.array(  [6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5 ])*1e-3, 
                            water_ec = 0.05)

water4 = Water(sample4WEC) #                  
print('water4', water4) #                      [0.05  0.11  0.08  0.11    nan   nan 0.067 0.07 ]
print('sample4WEC.Lw', sample4WEC.Lw)       

sample4WEC.info.to_excel('sample4WEC_info.xlsx')
sample4WEC.df.to_excel('sample4WEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example4')
#ax.plot(WunderlichC(sample4WEC.df.water, 0.006, 0.05, sample4WEC.df.water_ec, sample4WEC.Lw), sample4WEC.water, 'bo', alpha=0.3, markersize=8)
ax.plot(sample4WEC.bulk_ec, sample4WEC.water, 'bo',  markersize=4)
ax.plot(sample4WEC.bulk_ec, water4, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1*1e-3, 32*1e-3)
ax.set_ylim(-0.1, 0.65)
plt.show()


print("################## Example4b ####################") 
                         #                      0     1     2     3       4       5       6       7     8       9
sample4bWEC = Soil(water = np.array([           0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_ec=np.array([  6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_ec = 0.05
                            )

water4b = Water(sample4bWEC) #                
print('water4b', water4b) #                    
print('sample4bWEC.Lw', sample4bWEC.Lw) #      [0.05  0.11  0.08  0.11  0.079   nan 0.067 0.07  0.073 0.073]

sample4bWEC.info.to_excel('sample4bWEC_info.xlsx')
sample4bWEC.df.to_excel('sample4bWEC_df.xlsx')


##############
########################## Testing dielectric dispersion, fitting and non-fitting ###############################
##############

print("################## Example5 ####################")
                         #                     0     1     2     3     4     5     6      7    8     9      10
sample5WEC = Soil( bulk_ec=np.array([          10,   15,   20,   25,   7,    1,    12,    20,  5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_ec = 0.01,
                      frequency_ec = 500)

water5 = Water(sample5WEC) #                  
print("water5", water5) #                      
#                                              [0.566 0.65  0.65  0.65  0.421 0.15  0.639 0.65  0.316 0.65  0.65 ]
sample5WEC.info.to_excel('sample5WEC_info.xlsx')
sample5WEC.df.to_excel('sample5WEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example5')
ax.plot(sample5WEC.bulk_ec, [0.557, 0.65,  0.65,  0.65,  0.421, 0.074, 0.639, 0.65,  0.32,  0.65,  0.65 ], 'bo',  markersize=4)
ax.plot(sample5WEC.bulk_ec, [0.566, 0.65,  0.65,  0.65,  0.421, 0.15,  0.639, 0.65,  0.316, 0.65,  0.65 ], 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1*1e-3, 32*1e-3)
ax.set_ylim(-0.1, 0.65)
plt.show()


print("################## Example5b ####################")
                         #                    0      1     2       3         4      5     6      7      8     9      10
sample5bWEC = Soil( bulk_ec=np.array([           10,    0,    np.nan, np.nan,   7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.1,
                      frequency_ec = np.array([2e6]))

water5b = Water(sample5bWEC)  #             
print("water5b", water5b) #                   
#                                             [0.216 0.018   nan   nan 0.168 0.029 0.242 0.34  0.13  0.34  0.361]
sample5bWEC.info.to_excel('sample5bWEC_info.xlsx')
sample5bWEC.df.to_excel('sample5bWEC_df.xlsx')


print("################## Example6 ####################")     
                         #                      0     1     2     3     4     5     6     7     8     9     10
sample6WEC = Soil(water = np.array([            0.20, 0.30, 0.35]), 
                            bulk_ec=np.array([  10,   15,   20,   8.5,  8,    1,    12,   22,   5,    20,   30   ])*1e-3, 
                            bulk_density=1.7,
                            clay = 40,
                            water_ec=0.4,
                            frequency_ec=5e3
                            )

water6 = Water(sample6WEC) #                      [0.2   0.3   0.35  0.162 0.148 0.005 0.241 0.374 0.033 0.354 0.131]
print("water6", water6) #                         
print("sample6WEC.Lw", sample6WEC.Lw) #                         

sample6WEC.info.to_excel('sample6WEC_info.xlsx')
sample6WEC.df.to_excel('sample6WEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example6')
ax.plot(sample6WEC.bulk_ec, sample6WEC.water, 'bo',            markersize=4)
ax.plot(sample6WEC.bulk_ec, water6,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()


print("################## Example6b ####################")
                         #                     0     1      2      3      4      5     6      7      8     9      10
sample6bWEC = Soil( bulk_ec=np.array([         10,   15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                              water =         [0.1,  0.12],
                              bulk_density=1.7,
                              texture = 'Sand',
                              water_ec = 0.1,
                        frequency_ec=np.array([1,    2,     2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

water6b = Water(sample6bWEC) #               
print("water6b", water6b)               #        
                                              #  [0.1   0.12  0.394 0.447 0.213 0.065 0.293 0.394 0.173 0.394 0.406]
sample6bWEC.info.to_excel('sample6bWEC_info.xlsx')
sample6bWEC.df.to_excel('sample6bWEC_df.xlsx')
print('sample6bWEC', sample6bWEC.Lw)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example6b')
ax.plot(sample6bWEC.bulk_ec, sample6bWEC.water, 'bo',            markersize=4)
ax.plot(sample6bWEC.bulk_ec, water6b,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()


print("################## Example6c ####################")
                         #                    0      1       2     3      4      5     6      7      8     9      10
sample6cWEC = Soil( bulk_ec=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Clay',
                            water_ec = 0.1,
                    frequency_ec = np.array([ 1,    2,      2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

water6c = Water(sample6cWEC) #                  
print("water6c", water6c) #                     
#                                               [0.023 0.034 0.045 0.056 0.016 0.002 0.027 0.045 0.011 0.045 0.047]
sample6cWEC.info.to_excel('sample6cWEC_info.xlsx')
sample6cWEC.df.to_excel('sample6cWEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example6c')
ax.plot(sample6cWEC.bulk_ec, water6c, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()


print("################## Example7 ####################")
                        #                      0     1     2     3       4       5       6       7     8       9
sample7WEC = Soil(water =          [           0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                            bulk_ec =np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.05,
                    frequency_ec = [           50e1, 50e2, 50e2, 200e2,  200e2,  200e2,  50e2,   50e1, 50e1,   200e2] )

# bulk_ec_dc [0.0056  0.01014 0.00826     nan 0.00078     nan 0.00732 0.00798 0.00798  0.0076 ]

water7 = Water(sample7WEC)  #           
print('water7', water7) #                     [0.05  0.11  0.08  0.11    0.071   nan     0.028   0.07  0.028   0.028]

sample7WEC.info.to_excel('sample7WEC_info.xlsx')
sample7WEC.df.to_excel('sample7WEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example7')

ax.plot([0.0056, 0.01014, 0.00826, np.nan, 0.00078, np.nan, 0.00732, 0.00798, 0.00798,  0.0076], sample7WEC.df.water.values, '*',  markersize=4)

ax.plot(sample7WEC.df.bulk_ec, sample7WEC.water, 'bo',  markersize=4)
ax.plot(sample7WEC.df.bulk_ec, sample7WEC.df.water.values, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 20*1e-3)
ax.set_ylim(0, 0.30)
plt.show()


print("################## Example7b ####################")
                        #                       0     1     2     3       4       5       6       7     8       9
sample7bWEC = Soil(water = np.array([           0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                            bulk_ec=np.array([  6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.05,
                    frequency_ec =        [     50,   5,   50,    200e2,  200e2,  200e2,  50e2,   50e1, 50,     20]
                    )

water7b = Water(sample7bWEC)                  
print('water7b', water7b) #                    
                             #                 [0.05  0.11  0.08  0.11  0.071   nan 0.067 0.07  0.073 0.073]
sample7bWEC.info.to_excel('sample7bWEC_info.xlsx')
sample7bWEC.df.to_excel('sample7bWEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example7b')
ax.plot(sample7bWEC.df.bulk_ec, sample7bWEC.water, 'bo',  markersize=4)
ax.plot(sample7bWEC.df.bulk_ec, sample7bWEC.df.water.values, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32*1e-3)
ax.set_ylim(-0.1, 0.65)
plt.show()


print("################## ExampleSandValthe ####################")
                         #                     0     1     2     3     4      5      6      7     8     9      
samplevWEC = Soil( bulk_ec=np.array([          3,    8,    15,   20,   22,    7,     12,    18,   10,   2     ])*1e-3, 
                            bulk_density=1.4,
                            texture = 'Sand',
                            CEC = 1.6,
                            water_ec = 0.1,
                      frequency_ec = np.array([50,   500,  5000, 2000, 50000, 55000, 25000, 100,  10,   50]))

waterv = Water(samplevWEC) #                  [0.114 0.209 0.306 0.366 0.378  0.182  0.262  0.356 0.25  0.085]
print("waterv", waterv) #                        

samplevWEC.info.to_excel('samplevWEC_info.xlsx')
samplevWEC.df.to_excel('samplevWEC_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('ExampleSandValthe')
ax.plot(samplevWEC.bulk_ec, waterv, 'ro', alpha=0.3, markersize=8)

ax.set_xlim(0, 32*1e-3)
ax.set_ylim(0, 0.65)
plt.show()
