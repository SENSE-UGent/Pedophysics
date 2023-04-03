import matplotlib.pyplot as plt
import numpy as np

from pedophysics.predict import BulkCond, BulkPerm
from pedophysics.simulate import Soil

print("################## Example0 ####################")     
                         #                      0      1      2       3       4       5       6       7
sample0 = Soil( 
                            water_cond = 0.1,
                            frequency_cond=np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ]),
                            frequency_perm=np.array([0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ]))

print(np.all(sample0.df.frequency_perm == sample0.df.frequency_perm[0]))


print("################## Example1 ####################")     
                         #                      0      1      2       3       4       5       6       7
sample1 = Soil(water = np.array([   0.05,  0.1,   0.08,   0.11,   0.01,   np.nan, np.nan, 0.07    ]), 
                            bulk_cond=np.array([0.006, 0.011, 0.009,  np.nan, np.nan, np.nan, 0.008,  0.0085  ]), 
                            water_cond = 0.1,
                            instrument = 'TDR')

pred_bulk_cond1 = BulkCond(sample1) 
print("BulkCond(sample1) :", pred_bulk_cond1)
print("sample1.Lw :", sample1.Lw)
print(pred_bulk_cond1==np.array([0.006, 0.01127, 0.00921, 0.01228, np.nan, np.nan, 0.008, 0.00816]))

print("################### Example1b ###################")   
                         #                      0      1      2       3       4       5       6       7    
sample1b = Soil(water = np.array([  0.05,  0.1,   0.08,   0.11,   0.01,   np.nan, np.nan, 0.07    ]), 
                            bulk_cond=np.array([0.006, 0.011, 0.009,  np.nan, np.nan, np.nan, 0.008,  0.0085  ]), 
                            bulk_density=1.7,
                            water_cond = 0.1,
                            clay=2,
                       #     solid_perm = 5,
                       #     instrument= 'GPR')
                            frequency_cond = np.array([80, 80, 80, 80, 80, 80, 80, 80]))
pred_bulk_cond1b = BulkCond(sample1b) 
print("BulkCond(sample1b) :", pred_bulk_cond1b)
print("sample1b.Lw :", sample1b.Lw)
print(pred_bulk_cond1b==np.array([0.006, 0.01127, 0.00921, 0.01228, np.nan, np.nan, 0.008, 0.00816]))

print("################## Example1c ####################")     
# Silt Loam sample Wunderlich et al., 2013 #### 0      1      2       3       4       5       6       7       8
sample1c = Soil(water = np.array([  0.06,  0.08,  0.1,    0.12,   0.14,   np.nan, 0.17,   0.19,   0.28    ]), 
                            bulk_cond=np.array([2e-3, 5.5e-2, np.nan, np.nan, np.nan, np.nan, 1.2e-2, np.nan, np.nan  ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_cond = 0.1,
                            solid_perm = 5,
                            instrument = 'TDR')

pred_bulk_cond1c = BulkCond(sample1c) 
print("BulkCond(sample1c) :", pred_bulk_cond1c)
print("sample1c.Lw :", sample1c.Lw)
print(pred_bulk_cond1c==np.array([0.002, 0.00634, 0.00939, 0.01233, 0.01522, np.nan, 0.01953, 0.02239, 0.01102]))


print("############## Example2 #####################")       
                         #                      0     1    2     3       4       5       6       7
sample2 = Soil(water = np.array([   0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            instrument= 'GPR',
                frequency_perm = np.array([     2e9,  2e9, 2e9,  2e9,    2e9,    5e9]))
print("BulkPerm(sample2) :", BulkPerm(sample2)) 

print("############## Example2b #####################")       
                         #                      0     1    2     3       4       5       6       7
sample2b = Soil(water = np.array([  0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            instrument= 'GPR',
                            water_cond = 0.1,
                            clay = 10,
                            frequency_perm = np.array([2e9, 2e9, 2e9, 2e9, 2e9, 5e9]))

print("BulkPerm(sample2b) :", BulkPerm(sample2b)) 
print("sample2b.bulk_cond = ", sample2b.bulk_cond) 
print("sample2b.df.bulk_cond = ", sample2b.df.bulk_cond) 



print("############ Example3 ##################")
sample3 = Soil(water = 0.1, 
                            bulk_cond=np.array([ 0.0072,    0.007,   0.0075,  0.008]), 
                            sand=20.0,
                            silt = 40)
print("BulkCond(sample3) :", BulkCond(sample3))

print("############ Example3b ##################")
sample3b = Soil(water =              0.1, 
                            bulk_cond=np.array([ 0.0072,  0.007, 0.0075,  np.nan]), 
                            sand=20.0,
                            silt = 10,
                            bulk_density=1.5,
                            water_cond = 0.05,
                            instrument = 'GPR')

print("BulkCond(sample3b) :", BulkCond(sample3b))
print("sample3b.Lw :", sample3b.Lw)



print("############## Example4 ##############")
sample4 = Soil(water = np.array([   0.06,    0.08,   0.095,  0.128]), 
                            bulk_cond=np.array([0.007,   0.0072, 0.0075, np.nan]), 
                            temperature=25+273.15)

print("BulkCond(sample4) :", BulkCond(sample4))
print('***************np.nanmean(np.nan)', np.nanmean(np.nan))
print('***************np.nanmean([np.nan])', np.nanmean([np.nan]))

print("sample4.Lw :", sample4.Lw)

np.nanmean
# TODO there is a probleem here# aaaaahhhh there is not water_cond


print("################## Example5 ###################")
        #                                       0       1       2      3       4     5
sample5 = Soil(water = np.array([   0.06,   0.08,   0.095, 0.128             ]), 
                            bulk_cond=np.array([0.01,   0.014,  0.016, 0.02,   0.03, 0.04]), 
                            temperature=25.+273.15)

print("BulkCond(sample5) :", BulkCond(sample5))
# TODO que devuelvo cuando bulk_cond se pasa y la predigo, pero no hay cambios?



print(" ###############  Example6 #################")
sample6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25.+273.15,
                            clay = 20,
                            water_cond = 0.2, 
                            bulk_density=1.7,
                            frequency_cond = 60)
print("BulkCond(sample6) :", BulkCond(sample6))

print(" ###############  Example6b #################")
sample6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25.+273.15,
                            clay = 20,
                            water_cond = 0.2, 
                            bulk_density=1.7,
                            solid_cond=0.001)

print("BulkCond(sample6b) :", BulkCond(sample6b))



print("############### Example7 #################")
sample7 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25.+273.15, texture = "Clay", 
                            instrument = 'HydraProbe', 
                            pH_water=6, 
                            orgm = 0.4)

print('sample7.frequency_cond', sample7.frequency_cond)
print('sample7.df.frequency_cond', sample7.df.frequency_cond)
print("BulkCond(sample7) :", BulkCond(sample7))
print('sample7.df.frequency_cond', sample7.df.frequency_cond)
# TODO CEC module !!


print("############# Example8 ##################")
sample8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_density=1.6, 
                            sand = 20, 
                            silt = 60, 
                            water_cond = 0.05, 
                            CEC = 20., 
                            frequency_perm = 50e6)
print("BulkCond(sample8) :", BulkCond(sample8))
print('soil.Lw', sample8.Lw)



print("#############  Example9  #############")
sample9 = Soil(water = np.array([0.005,   0.01,   0.015, 0.02]), 
                            frequency_perm=np.array([50e6, 50e6, 50e6, 50e6]), 
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            pH_water = 7., 
                            CEC = 20.)
print("BulkCond(sample9) :", BulkCond(sample9)) 



print("#################  Example10  ####################")
sample10 = Soil(bulk_perm = np.array([np.nan, 7]),
                            frequency_cond=np.array([10 , 50,  100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 1e5]), 
                            water = np.array([       0.1, 0.1, 0.1, 0.1, 0.1, 0.1,  np.nan,  0.1,  0.1,   0.1,   0.1,   0.1]),
                            bulk_density=1.5,
                            water_cond = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print("BulkCond(sample10) :",BulkCond(sample10)) 


print("#################  Example11  ####################")
sample11 = Soil(bulk_perm = np.array([np.nan, 7]),
                            frequency_cond=np.array([1000, 1000, 1000, 1000, 1000, np.nan, 1000, 1000, 1000, 1000,   1000, 10]), 
                            water = np.array([       0.1,  0.1,  0.1,  0.1,  0.1,  0.1,    0.1,  0.1,  0.1,  np.nan, 0,    0.1]),
                            bulk_density=1.5,
                            water_cond = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print("BulkCond(sample11) :",BulkCond(sample11)) 


print("#################  Example12  ####################")
sample12 = Soil(bulk_perm = np.array([np.nan, 7]),
                            #                        0    1    2     3     4     5      6      7       8       9         10   11
                            frequency_cond=np.array([10,  20,  30,   55,   70,   66,    0,     np.nan, 48,     100,      5       ]), 
                            water = np.array([       0.1, 0.2, 0.11, 0.32, 0.61, 0.41,  0.01,  0.151,  0.21,   np.nan,   0,   0.1]),
                            bulk_density=1.5,
                            water_cond = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print("BulkCond(sample12) :",BulkCond(sample12)) # no hay bulk_density


print("#################  Example13  ####################")
sample13 = Soil(bulk_perm = np.array([np.nan, 7]),
                            frequency_cond=np.array([10,  10,  10,  10,  10,  np.nan, 10,  10,  10,  10,     10, 10]), 
                            water = np.array([       0.1, 0.1, 0.1, 0.1, 0.1, 0.1,    0.1, 0.1, 0.1, np.nan, 0,  0.1]),
                            bulk_density=1.5,
                            water_cond = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print("BulkCond(sample13) :",BulkCond(sample13)) # no hay bulk_density


print("############ Example14 ##################")
sample14 = Soil(water =   np.array([ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ]), 
                            bulk_cond=np.array([ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ]), 
                            frequency_cond=np.array([ 10]),
                            sand=20.0, silt = 10, bulk_density=1.5, water_cond = 0.05, instrument = 'GPR')

print("BulkCond(sample14) :", BulkCond(sample14))
print("sample14.Lw :", sample14.Lw)

print("############ Example14b ##################")
sample14b = Soil(water =   np.array([ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ]), 
                            bulk_cond=np.array([ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ]), 
                            frequency_cond=np.array([ 1000]),
                            sand=20.0, silt = 10, bulk_density=1.5, water_cond = 0.05, instrument = 'GPR')

print("BulkCond(sample14b) :", BulkCond(sample14b))
print("sample14b.Lw :", sample14b.Lw)

print("############ Example14c ##################")
sample14c = Soil(water =   np.array([ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ]), 
                            bulk_cond=np.array([ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014  ]), 
                            frequency_cond=np.array([ 10,  20,  30,   55,   70,   66,    0,     np.nan, 48,     99,      5       ]),
                            sand=20.0, silt = 10, bulk_density=1.5, water_cond = 0.05, instrument = 'GPR')

print("BulkCond(sample14c) :", BulkCond(sample14c))
print("sample14c.Lw :", sample14c.Lw)

print("############ Example14d ##################")
sample14d = Soil(water =   np.array([ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ]), 
                            bulk_cond=np.array([ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ]), 
                            frequency_cond=np.array([ 1000, 1500, 120, 150 ,200, 500, 101, 100, 800, 1200, 1e6, 1e8]),
                            sand=20.0, silt = 10, bulk_density=1.5, water_cond = 0.05, instrument = 'GPR')

print("BulkCond(sample14d) :", BulkCond(sample14d))
print("sample14d.Lw :", sample14d.Lw)

print("############ Example14e ##################")
sample14e = Soil(water =    np.array([ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ]), 
                            bulk_cond = np.array([ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ]), 
                            frequency_cond=np.array([ 1000]),
                            sand=20.0, silt = 10, bulk_density=1.5, water_cond = 0.05, instrument = 'GPR')

print("BulkCond(sample14e) :", BulkCond(sample14e))
print("sample14e.Lw :", sample14e.Lw)


print("#################  Graph example  ####################")
sample11g = Soil(water = np.array([0.1]), 
                            bulk_perm=np.array([7]), 
                            frequency_cond=np.logspace(2, 8, 100), 
                            temperature=25.+273.15, 
                            water_cond = 0.05,
                            bulk_density=1.5,
                            sand = 20, 
                            silt = 60)

sample12g = Soil(water = np.array([0.2]), 
                            bulk_perm=np.array([7]), 
                            frequency_cond=np.logspace(2, 8, 100), 
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            water_cond = 0.05,
                            bulk_density=1.5,
                            bulk_perm_inf = 5)

sample13g = Soil(water = np.array([0.3]), 
                            bulk_perm=np.array([7]), 
                            frequency_cond=np.logspace(2, 8, 100), 
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            water_cond = 0.05,
                            bulk_density=1.5)

sample14g = Soil(water = np.array([0.4]), 
                            bulk_perm=np.array([7]), 
                            frequency_cond=np.logspace(2, 8, 100), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            water_cond = 0.05,
                            bulk_density=1.5,
                            CEC = 20., 
                            bulk_perm_inf = 5)

#print('sample11.frequency_perm', sample11.frequency_perm)
#print('sample11.bulk_cond', sample11.bulk_cond)

pred_ls11g = BulkCond(sample11g)
pred_ls12g = BulkCond(sample12g)
pred_ls13g = BulkCond(sample13g)
pred_ls14g = BulkCond(sample14g)

plt.semilogx(sample11g.frequency_cond, pred_ls11g, 'bo', markersize=2)
plt.semilogx(sample12g.frequency_cond, pred_ls12g, 'ro', markersize=2)
plt.semilogx(sample13g.frequency_cond, pred_ls13g, 'go', markersize=2)
plt.semilogx(sample14g.frequency_cond, pred_ls14g, 'yo', markersize=4)

plt.show()