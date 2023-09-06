import numpy as np
import matplotlib
matplotlib.use('TkAgg',force=True)
from matplotlib import pyplot as plt
print("Switched to:",matplotlib.get_backend())

import sys
sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
#sys.path.insert(0, 'C:\\Users\\gasto\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import Water

##############
########################## Testing Fixed frequnecy, fitting and non-fitting ###############################
##############


print("################## Example0 ####################") 
                         #          0     1   2     3     4   5   6    7    8   9    10
sample0P = Soil( bulk_perm =        [10,   15, 20,   25,   7,  1,  12,  22,  5,  20,  30   ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'GPR')

water0 = Water(sample0P) # [0.129 0.21  0.283 0.35  0.074 0.    0.162 0.31  0.034 0.283 0.414]
print("water0", water0) # 

#sample0P.info.to_excel('sample0P_info.xlsx')
#sample0P.df.to_excel('sample0P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example0')
ax.plot(sample0P.bulk_perm, sample0P.water, 'bo',            markersize=4)
ax.plot(sample0P.bulk_perm, water0,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example0b ####################") 
                         #                      0     1     2     3       4       5       6       7     8       9
sample0bP = Soil(water =          [             0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                            bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ]), 
                            solid_perm = 5)

water0b = Water(sample0bP) #                    [0.05  0.11  0.08  0.11    nan     nan     nan     0.07  nan     nan]
print('water0b', water0b) # There is not a frequency_perm, then this is advised in the info.bulk_perm

#sample0bP.info.to_excel('sample0bP_info.xlsx')
#sample0bP.df.to_excel('sample0bP_df.xlsx')


print("################## Example1 ####################")     
                         #                       0     1     2     3       4       5       6       7
sample1P = Soil(water =                         [0.05, 0.11, 0.08, np.nan, np.nan, 0.07              ], 
                            bulk_perm=          [6,    11,   9,    np.nan, 8,      8.5,    12        ], 
                            bulk_density=1.7,
                            instrument = 'TDR')

water1 = Water(sample1P) #                      [0.05  0.11  0.08  nan     0.071   0.07    0.117     ]  
print('water1', water1) #                       

#sample1P.info.to_excel('sample1P_info.xlsx')
#sample1P.df.to_excel('sample1P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1')
ax.plot(sample1P.bulk_perm, sample1P.water, 'bo',  markersize=4)
ax.plot(sample1P.bulk_perm, water1, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1, 32)
ax.set_ylim(-0.1, 0.55)
plt.show()


print("################## Example1b ####################")     
                         #                      0     1     2     3       4       5       6       7       8       9
sample1bP = Soil(water =                       [0.05, 0.11, 0.08, np.nan, np.nan, 0.07,   np.nan, 0.2,    0.02,   np.nan ], 
                            bulk_perm=         [6,    11,   9,    np.nan, 8,      8.5,    14,     np.nan, np.nan, 1      ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'TDR')

water1b = Water(sample1bP) #                   [0.05  0.11  0.08  nan     0.071   0.07    0.194   0.2     0.02    0.   ]
print('water1b', water1b) #                    

#sample1bP.info.to_excel('sample1bP_info.xlsx')
#sample1bP.df.to_excel('sample1bP_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1b')
ax.plot(sample1bP.bulk_perm, sample1bP.water, 'bo',  markersize=4)
ax.plot(sample1bP.bulk_perm, water1b, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1, 32)
ax.set_ylim(-0.1, 0.55)
plt.show()


print("################## Example1c ####################")
      
                         #                      0     1     2     3     4    5  6   7   8  9   10
sample1cP = Soil(water =                       [0.20, 0.31, 0.36, 0.38, 0.05                        ], 
                            bulk_perm=         [10,   15,   20,   25,   7,   1, 12, 22, 5, 20, 30   ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5)

water1c = Water(sample1cP) #                   [0.2   0.31  0.36  0.38  0.05  nan  nan  nan  nan  nan  nan]
print("water1c", water1c) # 

#sample1cP.info.to_excel('sample1cP_info.xlsx')
#sample1cP.df.to_excel('sample1cP_df.xlsx')


print("################## Example3 ####################")  
                         #               0     1     2     3    4     5    6     7     8     9     10
sample3P = Soil(water = [                0.20, 0.30, 0.35                                             ], 
                            bulk_perm = [10,   15,   20,   8.5, 8,    1,   12,   22,   5,    20,   30   ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'GPR')

water3 = Water(sample3P) #              [0.2   0.3   0.35 0.173 0.164 0.   0.234 0.392 0.104 0.362 0.414]
print("water3", water3) #              

#sample3P.info.to_excel('sample3P_info.xlsx')
#sample3P.df.to_excel('sample3P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example3')
ax.plot(sample3P.bulk_perm, sample3P.water, 'bo',  markersize=4)
ax.plot(sample3P.bulk_perm, water3, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example4 ####################")     
                         #                      0      1     2    3      4    5    6     7     8     9     10
sample4P = Soil(water =                        [0.20,  0.30, 0.35                                               ], 
                            bulk_perm=         [10,    15,   20,  8.5,   8,   1,   12,   22,   5,    20,   30   ], 
                            bulk_density=1.7,
                            clay = 40,
                            solid_perm = 5,
                            instrument = 'TDR')

water4 = Water(sample4P)  #    [0.2   0.3   0.35  0.173  0.164 0.  0.234 0.392 0.104 0.362 0.459]
print("water4", water4) #      

#sample4P.info.to_excel('sample4P_info.xlsx')
#sample4P.df.to_excel('sample4P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example4')
ax.plot(sample4P.bulk_perm, sample4P.water, 'bo',  markersize=4)
ax.plot(sample4P.bulk_perm, water4, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example5 ####################")     
                         #               0     1     2   3    4  5  6   7   8  9   10
sample5P = Soil(water =                  [0.20, 0.30, 0.35                               ], 
                            bulk_perm = [10,   15,   20, 8.5, 8, 1, 12, 22, 5, 20, 30   ], 
                            bulk_density=1.7,
                            alpha = 0.3,
                            frequency_perm = [150e6])

water5 = Water(sample5P) # [0.2   0.3   0.35  0.173 0.164 0.    0.234 0.392 0.104 0.362 0.534]
print("water5", water5) # 

#sample5P.info.to_excel('sample5P_info.xlsx')
#sample5P.df.to_excel('sample5P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example5')
ax.plot(sample5P.bulk_perm, sample5P.water, 'bo',            markersize=4)
ax.plot(sample5P.bulk_perm, water5,        'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


##############
########################## Testing dielectric dispersion, fitting and non-fitting ###############################
##############


print("################## Example7 ####################")
                         #                      0     1     2     3       4       5       6       7     8       9
sample7P = Soil(             water =            [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                            bulk_perm =        [6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5   ], 
                            frequency_perm =   [50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6 ],
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.05)

water7 = Water(sample7P) #                      [0.05  0.11  0.08  0.11    0.05    nan     0.055   0.07  0.06     0.105]
print('water7', water7) # water7               [0.05  0.11  0.08  0.11    0.      nan     0.06    0.07  0.065    0.111]
#                                              [0.05  0.11  0.08  0.11    0.      nan     0.      0.07  0.    0.111] !!!
#                                              
#sample7P.info.to_excel('sample7P_info.xlsx')
#sample7P.df.to_excel('sample7P_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example7')
ax.plot(sample7P.df.bulk_ec, sample7P.water, 'bo',  markersize=4)
ax.plot(sample7P.df.bulk_ec, water7, 'ro', alpha=0.3, markersize=8)

ax.plot(sample7P.bulk_perm, sample7P.water, 'bo',  markersize=4)
ax.plot(sample7P.bulk_perm, water7, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)

plt.show()


print("################## Example7b ####################") 
       
                         #                      0     1     2     3       4       5       6       7
sample7bP = Soil(water =                       [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                            bulk_perm =        [6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5   ], 
                            frequency_perm =   [50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6 ],
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5)

water7b = Water(sample7bP)  #                  
print('water7b', water7b) #                    [0.05 0.11 0.08 0.11  nan  nan  nan 0.07  nan  nan] no hay water ec

#sample7bP.info.to_excel('sample7bP_info.xlsx')
#sample7bP.df.to_excel('sample7bP_df.xlsx')


print("################## Example7c ####################") 
       
                         #      0       1       2       3       4       5       6       7       8       9
sample7cP = Soil(water =       [0.05,   0.11,   0.08,   0.11,   np.nan, np.nan, np.nan, 0.07,   np.nan, np.nan], 
            bulk_perm =        [6,      11,     9,      np.nan, 1,      np.nan, 8,      8.5,    8.5,    8.5   ],
            bulk_ec =          [np.nan, np.nan, np.nan, 0.002,  np.nan, np.nan, 0.003, np.nan, np.nan, np.nan],
            frequency_perm =   [50e6,   50e6,   50e6,   200e6,  200e6,  200e6,  50e6,   50e6,   50e6,   200e6 ],
            frequency_ec =     [50e2,   50e2,   50e2,   200e1,  200e3,  200e4,  50e3,   50e3,   50e3,   20    ],
                              #[0.05    0.11    0.08    0.11    0.05    nan     0.055   0.07    0.06    0.105 ]
                              #[0.05    0.11    0.08    0.11    0.15    nan     0.15    0.07    0.15    0.15  ]
                              #[0.05    0.11    0.08    0.11    0.      nan     0.06    0.07    0.065   0.111 ]
                              #[0.05    0.11    0.08    0.11    0.15    nan     0.15    0.07    0.15    0.15  ]
                              #[0.05    0.11    0.08    0.11    0.      nan     0.111   0.07    0.      0.117]
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.05)

water7c = Water(sample7cP)  
print('water7c', water7c) 

#sample7cP.info.to_excel('sample7cP_info.xlsx')
#sample7cP.df.to_excel('sample7cP_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example7c')
ax.plot(sample7cP.bulk_perm, sample7cP.water, 'bo',  markersize=4)
ax.plot(sample7cP.bulk_perm, water7c, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(-1, 32)
ax.set_ylim(-0.1, 0.55)
plt.show()


print("################## Example8 ####################")
                         #                       0     1      2      3      4      5     6      7      8       9      10
sample8P = Soil( bulk_perm =                    [10,   15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        solid_perm = 5,
                        frequency_perm =        [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6])

water8 = Water(sample8P) # [nan nan nan nan nan nan nan nan nan nan nan]
print("water8", water8)

#sample8P.info.to_excel('sample8P_info.xlsx')
#sample8P.df.to_excel('sample8P_df.xlsx')


print("################## Example8b ####################")
                         #                       0      1       2     3      4      5     6      7      8       9      10
sample8bP = Soil( bulk_perm =                   [10,    15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.1,
                      frequency_perm =         [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6])
#                                              
water8b = Water(sample8bP)  #                  [0.002 0.027  0.079  0.167  0.001  0.    0.1    0.437  nan     0.65   0.65 ]
print("water8b", water8b) #                    [0.019 0.034  0.081  0.169  0.018  0.018 0.105  0.437  nan     0.65   0.65 ]

#sample8bP.info.to_excel('sample8bP_info.xlsx')
#sample8bP.df.to_excel('sample8bP_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example8b')
ax.plot(sample8bP.bulk_perm, water8b, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)
plt.show()


print("################## Example8c ####################")
                         #                      0      1       2     3      4      5     6      7        8       9      10
sample8cP = Soil( bulk_perm =                  [10,    15,    20,    25,    7,     1,    12,    np.nan,  5,      20,    22 ], 
                            bulk_density=1.7,
                            texture = 'Clay',
                            solid_perm = 5,
                            water_ec = 0.1,
                        frequency_perm =       [1e6,   2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,    np.nan, 100e6, 200e6])

water8c = Water(sample8cP)                    # [0.    0.001 0.003 0.011  .    0.    0.005  nan   nan 0.272 0.649]
print("water8c", water8c) #                     [0.    0.001 0.003 0.011 0.    0.    0.005  nan   nan 0.272 0.65 ]

#sample8cP.info.to_excel('sample8cP_info.xlsx')
#sample8cP.df.to_excel('sample8cP_df.xlsx')

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example8c')
ax.plot(sample8cP.bulk_perm, water8c, 'ro', alpha=0.3, markersize=8)
ax.set_xlim(0, 32)
ax.set_ylim(0, 0.55)
plt.show()


print("################## Example9 ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample9P = Soil( bulk_perm =                    [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            water_ec = 0.1,
                            frequency_perm = 20e6)

water9 = Water(sample9P)  #                     [0.041 0.197 0.376 0.623 0.003 0.    0.083 0.376 0.    0.376 0.457]
print("water9", water9) #                       [0.044 0.196 0.376 0.623 0.019 0.018 0.085 0.376 0.018 0.376 0.457]

#sample9P.info.to_excel('sample9P_info.xlsx')
#sample9P.df.to_excel('sample9P_df.xlsx')


print("################## Example9b ####################")
                         #                      0      1       2     3      4      5     6      7      8     9      10
sample9bP = Soil( bulk_perm =                  [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
                              bulk_density=1.7,
                              texture = 'Sand',
                              solid_perm = 5,
                              water_ec = 0.1,
                              frequency_perm = 1e6)

water9b = Water(sample9bP) #                   [0.002 0.013 0.034 0.069 0.    0.    0.003 0.034 0.    0.034 0.045]
print("water9b", water9b) #                    [0.019 0.025 0.039 0.073 0.018 0.018 0.019 0.039 0.018 0.039 0.055]

#sample9bP.info.to_excel('sample9bP_info.xlsx')
#sample9bP.df.to_excel('sample9bP_df.xlsx')


print("################## ExampleSandValthe ####################")
                         #                      0     1        2      3      4      5     6      7      
samplevP = Soil( bulk_perm =                   [3,    8,       15,    20,    22,    7,    12,    18     ], 
                            bulk_density=1.4,
                            texture = 'Sand',
                            solid_perm = 5,
                            CEC = 1.6,
                      frequency_perm = 50e6  )

waterv = Water(samplevP) #                     [0.005 0.148 0.282 0.359 0.386 0.124 0.23  0.329]
print("waterv", waterv) #                   

#samplevP.info.to_excel('samplevP_info.xlsx')
#samplevP.df.to_excel('samplevP_df.xlsx')