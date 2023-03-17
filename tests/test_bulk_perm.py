from pedophysics.pedophysical_models.bulk_perm import LongmireSmithP
from pedophysics.predict import BulkPerm
from pedophysics.simulate import Soil

import numpy as np
import matplotlib.pyplot as plt


print("################## Example1 ####################")     
                         #                      0     1    2     3       4       5       6       7
sample1 = Soil(water = np.array([   0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            texture = 'Sand',
                            solid_perm = 5,
                            instrument = 'TDR')

print("predict.BulkPerm(sample1) :", BulkPerm(sample1))
print("predict.BulkPerm(sample1) :", BulkPerm(sample1))



print("################### Example1b ###################")   
                         #                      0     1    2     3       4       5       6       7
sample1b = Soil(water = np.array([  0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            clay=2,
                            solid_perm = 5,
                            instrument= 'GPR')

print("BulkPerm(sample1b) :", BulkPerm(sample1b))



print("############## Example2 #####################")       
                         #                      0     1    2     3       4       5       6       7
sample2 = Soil(water = np.array([   0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            instrument= 'GPR',
                            frequency_perm = np.array([2e9, 2e9, 2e9, 2e9, 2e9, 5e9]))

print('sample2.frequency_perm', sample2.frequency_perm)
print("BulkPerm(sample2) :", BulkPerm(sample2))   # Frequency changing module cannot work cuse bulk_cond is not finished



print("############ Example3 ##################")
sample3 = Soil(water = 0.1, 
                            bulk_perm=np.array([ 7.2,    7,   7.5,     8]), 
                            sand=20.0)

print("BulkPerm(sample3) :", BulkPerm(sample3))



print("############ Example3b ##################")
sample3b = Soil(water = 0.1, 
                            bulk_perm=np.array([ 7.2,    7,   7.5,     np.nan]), 
                            sand=20.0,
                            silt = 10,
                            instrument = 'GPR')

print("BulkPerm(sample3b) :", BulkPerm(sample3b))



print("############## Example4 ##############")
sample4 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_perm=np.array([7,      7.2,    7.5 ,  np.nan]), 
                            temperature=25.)
print("BulkPerm(sample4) :", BulkPerm(sample4))



print("################## Example5 ###################")
        #                                          0       1        2      3  4   5
sample5 = Soil(water = np.array([   0.06,   0.08,   0.095, 0.128       ]), 
                            bulk_perm=np.array([   7,    7.2,     7.5,     8, 9, 10]), 
                            temperature=25.)

print("predict.BulkPerm(sample5) :", BulkPerm(sample5))
# TODO que devuelvo cuando bulk_perm se pasa y la predigo, pero no hay cambios?



print(" ###############  Example6 #################")
sample6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25., 
                            frequency_perm = 60e6)
print("BulkPerm(sample6) :", BulkPerm(sample6))



print("############### Example7 #################")
sample7 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25., texture = "Clay", 
                            instrument = 'HydraProbe', 
                            pH_water=6, 
                            orgm = 0.4)

print('sample7.frequency_perm', sample7.frequency_perm)
print("BulkPerm(sample7) :", BulkPerm(sample7))
#print('sample7', sample7)
# TODO CEC module !!


print("############# Example8 ##################")
sample8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_density=1.6, 
                            sand = 20, 
                            silt = 60, 
                            pH_water = 7., 
                            CEC = 20., 
                            frequency_perm = 50e6)
print("BulkPerm(sample8) :", BulkPerm(sample8))
print('soil.Lw', sample8.Lw)



print("#############  Example9  #############")
sample9 = Soil(water = np.array([0.005,   0.01,   0.015, 0.02]), 
                            frequency_perm=np.array([50*10**6, 50*10**6, 50*10**6, 50*10**6]), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            pH_water = 7., 
                            CEC = 20.)
print("BulkPerm(sample9) :",BulkPerm(sample9)) # no hay bulk_density



print("#################  Example10  ####################")
sample10 = Soil(bulk_cond = np.array([0.0128, 0.0128, 0.0128, 0.0128, 0.0128, 0.0128]), 
                            bulk_perm = np.array([np.nan, 7]),
                            frequency_perm=np.array([10**6 , 10*10**6, 50*10**6, 100*10**6, 500*10**6, 1*10**9]), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
print('pdm.LongmireSmithP(sample10.bulk_cond, sample10.bulk_perm_inf, sample10.frequency_perm)', LongmireSmithP(sample10.bulk_cond, sample10.bulk_perm_inf, sample10.frequency_perm))
print("BulkPerm(sample10) :", BulkPerm(sample10))



print("#################  Graph example  ####################")
sample11 = Soil(bulk_cond = np.array([0.001]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample12 = Soil(bulk_cond = np.array([0.010]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample13 = Soil(bulk_cond = np.array([0.050]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample14 = Soil(bulk_cond = np.array([0.080]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(6, 8, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

print('sample11.frequency_perm', sample11.frequency_perm)
print('sample11.bulk_cond', sample11.bulk_cond)

pred_ls11 = LongmireSmithP(sample11.bulk_cond, sample11.bulk_perm_inf, sample11.frequency_perm)
pred_ls12 = LongmireSmithP(sample12.bulk_cond, sample13.bulk_perm_inf, sample12.frequency_perm)
pred_ls13 = LongmireSmithP(sample13.bulk_cond, sample13.bulk_perm_inf, sample13.frequency_perm)
pred_ls14 = LongmireSmithP(sample14.bulk_cond, sample14.bulk_perm_inf, sample14.frequency_perm)

plt.semilogx(sample11.frequency_perm, pred_ls11, 'bo', markersize=2)
plt.semilogx(sample12.frequency_perm, pred_ls12, 'ro', markersize=2)
plt.semilogx(sample13.frequency_perm, pred_ls13, 'go', markersize=2)
plt.semilogx(sample14.frequency_perm, pred_ls14, 'yo', markersize=4)

plt.show()