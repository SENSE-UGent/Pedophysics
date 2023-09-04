import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.pedophysical_models.bulk_perm import LongmireSmithP
from pedophysics.predict import BulkPerm
from pedophysics.simulate import Soil

##############
########################## Testing fixed frequency, fitting and non-fitting ###############################
##############


print("################## Example0 ####################")     
                         #   0     1    2     3       4       5       6       7
sampleP0 = Soil( water =     [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07], 
                bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                frequency_perm = [7e8, 7e8],
                bulk_density=1.7,
                solid_perm = 5,
                instrument = 'HydraProbe',
                salinity = 0.1)

print("predict.BulkPerm(sampleP0) :", BulkPerm(sampleP0)) 
print("predict.BulkPerm(sampleP0) :", BulkPerm(sampleP0))

sampleP0.info.to_excel('sampleP0_info.xlsx')
sampleP0.df.to_excel('sampleP0_df.xlsx')


print("################## Example1 ####################")     
                         #   0     1    2     3       4       5       6       7
sampleP1 = Soil( water =     [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07], 
                bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                bulk_density=1.7,
                solid_perm = 5,
                instrument = 'TDR')

print("predict.BulkPerm(sampleP1) :", BulkPerm(sampleP1)) 
print("predict.BulkPerm(sampleP1) :", BulkPerm(sampleP1))

sampleP1.info.to_excel('sampleP1_info.xlsx') # Fittea pero solo el estado 4 puede ser predicho, pq los demas estan fuera del rango aunque los intenta predecir con CRIM
sampleP1.df.to_excel('sampleP1_df.xlsx')


print("################### Example1b ###################")   
                         #               0     1    2     3       4       5       6       7
sampleP1b = Soil(water =                [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07]    , 
                            bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                            bulk_density=1.7,
                            clay=2,
                            solid_perm = 5,
                            instrument= 'GPR')
#                                       [6.    11.  9.    10.653  4.648   nan     8.      8.5  ]
print("BulkPerm(sampleP1b) :", BulkPerm(sampleP1b))

sampleP1b.info.to_excel('sampleP1b_info.xlsx')
sampleP1b.df.to_excel('sampleP1b_df.xlsx')


print("############ Example3 ##################")

sampleP3 = Soil( water = 0.1, 
                bulk_perm = [7.2,    7,   7.5,     8], 
                sand=20.0)
                          # [7.2     7.   7.5      8. ]
print("BulkPerm(sampleP3) :", BulkPerm(sampleP3))

sampleP3.info.to_excel('sampleP3_info.xlsx')
sampleP3.df.to_excel('sampleP3_df.xlsx')


print("############ Example3b ##################")
sampleP3b = Soil(water = 0.1, 
                bulk_perm =       [7.2,  7,   7.5,   np.nan], 
                sand = 20.0,
                silt = 10,
                instrument = 'GPR')
                                # [7.2   7.   7.5    7. ]
print("BulkPerm(sampleP3b) :", BulkPerm(sampleP3b))

sampleP3b.info.to_excel('sampleP3b_info.xlsx')
sampleP3b.df.to_excel('sampleP3b_df.xlsx')


print("############## Example4 ##############")
sampleP4 = Soil( water =     [0.06,   0.08,   0.095, 0.128], 
                bulk_perm = [7,      7.2,    7.5 ,  np.nan], 
                temperature=25.)
                          # [7.      7.2     7.5    nan]
print("BulkPerm(sampleP4) :", BulkPerm(sampleP4))

sampleP4.info.to_excel('sampleP4_info.xlsx')
sampleP4.df.to_excel('sampleP4_df.xlsx')


print(" ###############  Example6 #################")
sampleP6 = Soil(water = [0.06,   0.08,   0.095, 0.128], 
                temperature=25., 
                frequency_perm = 60e6)
# [nan nan nan nan]
print("BulkPerm(sampleP6) :", BulkPerm(sampleP6))

sampleP6.info.to_excel('sampleP6_info.xlsx')
sampleP6.df.to_excel('sampleP6_df.xlsx')


print("############### Example6b #################")
sampleP6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25., texture = "Clay", 
                            instrument = 'HydraProbe', CEC = 2, bulk_density = 1.3,
                            orgm = 0.4)
 
print('sampleP6b.frequency_perm', sampleP6b.frequency_perm)
print("BulkPerm(sampleP6b) :", BulkPerm(sampleP6b)) # [4.115 4.817 5.377 6.707]

sampleP6b.info.to_excel('sampleP6b_info.xlsx')
sampleP6b.df.to_excel('sampleP6b_df.xlsx')


print("############# Example6c ##################")
sampleP6c = Soil(water =  [0.06,  0.08,   0.095,  0.128], 
        bulk_density=1.6, 
        sand = 20, 
        silt = 60, 
        CEC = 20., 
        frequency_perm = 50e6)

print("BulkPerm(sampleP6c) :", BulkPerm(sampleP6c)) # [ 8.661 10.462 11.791 14.661]
print('soil.Lw', sampleP6c.Lw)

sampleP6c.info.to_excel('sampleP6c_info.xlsx')
sampleP6c.df.to_excel('sampleP6c_df.xlsx')


print("################## ExampleSandValthe ####################")
                        # 0        1       2       3        4        5        6        7      
samplePv = Soil( water =  [0.03,    0.08,   0.15,   0.20,    0.22,    0.07,    0.12,    0.18     ] , 
                bulk_density=1.4,
                texture = 'Sand',
                solid_perm = 5,
                CEC = 1.6,
                frequency_perm = np.array([50e6]))
                        # [3.697   5.315   8.092   10.447   11.476   4.967    6.828    9.468]

print("BulkPerm(samplePv)", BulkPerm(samplePv))

samplePv.info.to_excel('samplePv_info.xlsx')
samplePv.df.to_excel('samplePv_df.xlsx')


print("#################  Graph example  ####################")
sample11 = Soil(bulk_ec = np.array([0.001]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample12 = Soil(bulk_ec = np.array([0.010]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample13 = Soil(bulk_ec = np.array([0.050]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample14 = Soil(bulk_ec = np.array([0.080]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

pred_ls11 = LongmireSmithP(sample11.bulk_ec, sample11.bulk_perm_inf, sample11.frequency_perm)
pred_ls12 = LongmireSmithP(sample12.bulk_ec, sample13.bulk_perm_inf, sample12.frequency_perm)
pred_ls13 = LongmireSmithP(sample13.bulk_ec, sample13.bulk_perm_inf, sample13.frequency_perm)
pred_ls14 = LongmireSmithP(sample14.bulk_ec, sample14.bulk_perm_inf, sample14.frequency_perm)

plt.semilogx(sample11.frequency_perm, pred_ls11, 'bo', markersize=2)
plt.semilogx(sample12.frequency_perm, pred_ls12, 'ro', markersize=2)
plt.semilogx(sample13.frequency_perm, pred_ls13, 'go', markersize=2)
plt.semilogx(sample14.frequency_perm, pred_ls14, 'yo', markersize=4)

plt.show()


print("############## Example7 #####################")       
                         #        0     1    2     3       4       5       6       7
sampleP7 = Soil(water =          [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07] , 
                bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ] , 
                bulk_density=1.7,
                instrument= 'GPR',
                frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9]  )
#                                [6.    11.  9.    nan     nan     nan  8.   8.5]
print("BulkPerm(sampleP7) :", BulkPerm(sampleP7))
print("sampleP7.df.water_ec", sampleP7.df.water_ec)

sampleP7.info.to_excel('sampleP7_info.xlsx')
sampleP7.df.to_excel('sampleP7_df.xlsx')


print("############## Example7b #####################")       
                         #         0     1    2     3       4       5       6       7
sampleP7b = Soil(water =          [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07] , 
                bulk_perm =       [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ] , 
                bulk_density=1.7,
                instrument= 'GPR',
                water_ec = 0.2,
                clay = 10,
                frequency_perm =  [2e9,  2e9, 2e9,  2e9,    2e9,    5e9]  )
#                                 [6.    11.  9.    6.785   5.657   nan     8.      8.5  ]
print("BulkPerm(sampleP7b) :", BulkPerm(sampleP7b))

sampleP7b.info.to_excel('sampleP7b_info.xlsx')
sampleP7b.df.to_excel('sampleP7b_df.xlsx')


print("#################  Example8  ####################")
sampleP8 = Soil(bulk_ec =                     [0.0128, 0.0128, 0.0128, 0.0128, 0.0128, 0.0128], 
                            bulk_perm =      [np.nan, 7] ,
                            frequency_perm = [1e6 ,   1e6,    50e6,   100e6,  500e6,  1e9], 
                            #                [47.582  7.      14.049  12.271  8.888   8.134]
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print('pdm.LongmireSmithP', LongmireSmithP(sampleP8.bulk_ec, sampleP8.bulk_perm_inf, sampleP8.frequency_perm))
print("BulkPerm(sampleP8) :", BulkPerm(sampleP8))

sampleP8.info.to_excel('sampleP8_info.xlsx')
sampleP8.df.to_excel('sampleP8_df.xlsx')