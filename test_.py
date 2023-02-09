
import pedophysics
import numpy as np
from predict.bulk_perm import bulk_perm
import pedophysical_models as pdm
import matplotlib.pyplot as plt


print("Example1")        #                      0     1    2     3       4       5       6       7
sample1 = pedophysics.Soil(water = np.array([   0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            instrument= 'GPR')

print('sample1.frequency_perm', sample1.frequency_perm)
print("predict.bulk_perm(sample1) :", bulk_perm(sample1))


print("Example2")        #                      0     1    2     3       4       5       6       7
sample2 = pedophysics.Soil(water = np.array([   0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07]), 
                            bulk_perm=np.array([6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ]), 
                            bulk_density=1.7,
                            instrument= 'GPR',
                            frequency_perm = np.array([2e9, 2e9, 2e9, 2e9, 2e9, 5e9]))

print('sample2.frequency_perm', sample2.frequency_perm)
print("predict.bulk_perm(sample2) :", bulk_perm(sample2))



print("Example3")
sample3 = pedophysics.Soil(water = np.array([0.08, 0.06, 0.095, 0.128, 0.1, 0.2]), 
                            bulk_perm=np.array([7.2, 7, 7.5, 8]), 
                            sand=20.0)

#print("predict.bulk_perm(sample3) :", predict.bulk_perm(sample3))


print("Example4")
sample4 = pedophysics.Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_perm=np.array([7,      7.2,    7.5 ,  8]), 
                            temperature=25.)
#print("predict.bulk_perm(sample4) :", predict.bulk_perm(sample4))


print("Example5")
sample5 = pedophysics.Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_perm=np.array([7,      7.2,    7.5 ,  8, 9, 10]), 
                            temperature=25.)
print(len(sample5.frequency_perm))
#print("predict.bulk_perm(sample5) :", predict.bulk_perm(sample5))


print("Example6")
sample6 = pedophysics.Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25., 
                            frequency_perm = 60e6)
#print("predict.bulk_perm(sample6) :", predict.bulk_perm(sample6))


print("Example7")
sample7 = pedophysics.Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_perm=np.array([7]), 
                            temperature=25., texture = "Clay", 
                            instrument = 'HydraProbe', 
                            pH_water=6, 
                            orgm = 0.4)
#print('sample7.frequency_perm', sample7.frequency_perm)
#print("predict.bulk_perm(sample7) :", predict.bulk_perm(sample7))
#print('sample7', sample7)


print("Example8")
sample8 = pedophysics.Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_perm=np.array([7]), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            pH_water = 7., 
                            CEC = 20., 
                            frequency_perm = 50e6)
#print("predict.bulk_perm(sample8) :", predict.bulk_perm(sample8))
#print(sample8)


print("Example9")
sample9 = pedophysics.Soil(water = np.array([0.005,   0.01,   0.015, 0.02]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.array([50*10**6, 50*10**6, 50*10**6, 50*10**6]), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            pH_water = 7., 
                            CEC = 20.)
#print('bulk_perm_inf == []', (sample9.bulk_perm_inf))
#print('(all(soil.frequency_perm) <= 3e10)', (all(sample9.frequency_perm) <= 3e10))
#print('len(soil.frequency_perm) == len(soil.bulk_cond)', len(sample9.frequency_perm) == len(sample9.water))

#print("predict.bulk_perm(sample9) :", predict.bulk_perm(sample9))
#print(sample9)


print("Example10")
sample10 = pedophysics.Soil(bulk_cond = np.array([0.0128, 0.0128, 0.0128, 0.0128, 0.0128, 0.0128]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.array([10**6 , 10*10**6, 50*10**6, 100*10**6, 500*10**6, 1*10**9]), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
#print(pdm.longmire_smith(sample10.bulk_cond, sample10.bulk_perm_inf, sample10.frequency_perm))
#print("predict.bulk_perm(sample10) :", predict.bulk_perm(sample10))

sample11 = pedophysics.Soil(bulk_cond = np.array([0.001]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(7, 10, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample12 = pedophysics.Soil(bulk_cond = np.array([0.010]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(7, 10, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample13 = pedophysics.Soil(bulk_cond = np.array([0.030]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(7, 10, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

sample14 = pedophysics.Soil(bulk_cond = np.array([0.060]), 
                            bulk_perm=np.array([7]), 
                            frequency_perm=np.logspace(7, 10, 100)*2, 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

print('sample11.frequency_perm', sample11.frequency_perm)
print('sample11.bulk_cond', sample11.bulk_cond)

pred_ls11 = pdm.longmire_smith(sample11.bulk_cond, sample11.bulk_perm_inf, sample11.frequency_perm)
pred_ls12 = pdm.longmire_smith(sample12.bulk_cond, sample13.bulk_perm_inf, sample12.frequency_perm)
pred_ls13 = pdm.longmire_smith(sample13.bulk_cond, sample13.bulk_perm_inf, sample13.frequency_perm)
pred_ls14 = pdm.longmire_smith(sample14.bulk_cond, sample14.bulk_perm_inf, sample14.frequency_perm)

plt.semilogx(sample11.frequency_perm, pred_ls11, 'bo', markersize=2)
plt.semilogx(sample12.frequency_perm, pred_ls12, 'ro', markersize=2)
plt.semilogx(sample13.frequency_perm, pred_ls13, 'go', markersize=2)
plt.semilogx(sample14.frequency_perm, pred_ls14, 'yo', markersize=4)

plt.show()