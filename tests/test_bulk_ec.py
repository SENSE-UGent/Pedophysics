import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.predict import BulkEC, BulkPerm
from pedophysics.simulate import Soil

from pedophysics.pedophysical_models.bulk_ec import Rhoades


##############
########################## Testing DC frequency, fitting and non-fitting ###############################
##############

def test_sampleC0():
        sampleC0 = Soil(water = 0.1, 
                        bulk_ec= [ 0.0072,    0.007,   0.0075,  0.008], 
                        sand=20.0,
                        silt = 40)

        assert BulkEC(sampleC0) == np.array([0.0072, 0.007,  0.0075, 0.008])


#sampleC0.info.to_excel('sampleC0_info.xlsx')
#sampleC0.df.to_excel('sampleC0_df.xlsx')


print("############ Example0b ##################")
sampleC0b = Soil(water =              0.1, 
                bulk_ec = [ 0.0072,  0.007, 0.0075,  np.nan], 
                sand=20.0,
                silt = 10,
                bulk_density=1.5,
                water_ec = 0.05,
                instrument = 'GPR')
                                   #  [0.0072    0.007  0.0075   0.007 ]
print("BulkEC(sampleC0b) :", BulkEC(sampleC0b))
print("sampleC0b.Lw :", sampleC0b.Lw)
#sampleC0b.info.to_excel('sampleC0b_info.xlsx')
#sampleC0b.df.to_excel('sampleC0b_df.xlsx')


print("#################  Example0c  ####################")
sampleC0c = Soil(bulk_perm =                [np.nan, 7],
                frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                bulk_density=1.5,
                water_ec = 0.05,
                sand = 20, 
                silt = 60, 
                CEC = 20)

# [0.008388 0.008388 0.008815 0.008867 0.008924 0.008988      nan 0.008388 0.009239 0.009355 0.009528 0.009774]
print("BulkEC(sampleC0c) :",BulkEC(sampleC0c)) 
#sampleC0c.info.to_excel('sampleC0c_info.xlsx')
#sampleC0c.df.to_excel('sampleC0c_df.xlsx')


print("#################  Example0d  ####################")
sampleC0d = Soil(bulk_perm =                [np.nan, 7],
                frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                bulk_density=1.5,
                water_ec = 0.05,
                texture = 'Silt loam',
                instrument = 'EMI Dualem')

# [0.006302 0.006302 0.006654 0.006691 0.006739 0.006795      nan 0.006991 0.007006 0.007094 0.007257 0.007474]
print("BulkEC(sampleC0d) :",BulkEC(sampleC0d)) 
#sampleC0d.info.to_excel('sampleC0d_info.xlsx')
#sampleC0d.df.to_excel('sampleC0d_df.xlsx')


print("################## Example1 ####################")     
                         #                    0      1        2         3         4       5       6       7
sampleC1 = Soil(water =                      [0.05,  0.1,     0.08,     0.11,     0.01,   np.nan, np.nan, 0.07    ], 
                            bulk_ec=         [0.006, 0.011,   0.009,    np.nan,   np.nan, np.nan, 0.008,  0.0085  ], 
                            water_ec = 0.1,
                            instrument = 'TDR')

pred_bulk_ec1 = BulkEC(sampleC1)  #           [0.006  0.011    0.009     0.012277  nan     nan     0.008   0.0085  ]
print("BulkEC(sampleC1) :", pred_bulk_ec1)
print("sampleC1.Lw :", sampleC1.Lw)
print("sampleC1.df.water_ec :", sampleC1.df.water_ec)
#sampleC1.info.to_excel('sampleC1_info.xlsx')
#sampleC1.df.to_excel('sampleC1_df.xlsx')

print("################### Example1b ###################")   
                         #                      0      1        2        3        4       5       6       7    
sampleC1b = Soil(water =                        [0.05,  0.1,     0.08,    0.11,    0.01,   np.nan, np.nan, 0.07    ], 
                            bulk_ec =          [0.006, 0.011,   0.009,   np.nan,  np.nan, np.nan, 0.008,  0.0085  ], 
                            bulk_density=1.7,
                            water_ec = 0.1,
                            clay=2,
                            frequency_ec = 80)

pred_bulk_ec1b = BulkEC(sampleC1b) 
print("BulkEC(sample1b) :", pred_bulk_ec1b)#   [0.006  0.011    0.009    0.012277 0.000123 nan     0.008    0.0085  ]
print("sample1b.Lw :", sampleC1b.Lw)     
#sampleC1b.info.to_excel('sampleC1b_info.xlsx')
#sampleC1b.df.to_excel('sampleC1b_df.xlsx')

print("################## Example1c ####################")     
# Silt Loam sample Wunderlich et al., 2013 #### 0      1       2        3        4       5       6       7        8
sampleC1c = Soil(water =                      [  0.06,  0.08,   0.1,     0.12,    0.14,   np.nan, 0.17,   0.19,    0.28    ], 
                            bulk_ec=         [  2e-3,  5.5e-2, np.nan,  np.nan,  np.nan, np.nan, 1.2e-2, np.nan,  np.nan  ], 
                            bulk_density=1.7,
                            texture = 'Sand',
                            water_ec = 0.1,
                            solid_perm = 5,
                            instrument = 'TDR')

pred_bulk_ec1c = BulkEC(sampleC1c)  
                                           # [  0.002  0.055   0.009393 0.012327 0.01522  nan     0.012   0.022389 0.011016]
print("BulkEC(sample1c) :", pred_bulk_ec1c)
print("sample1c.Lw :", sampleC1c.Lw)
#sampleC1c.info.to_excel('sampleC1c_info.xlsx')
#sampleC1c.df.to_excel('sampleC1c_df.xlsx')


print("############## Example4 ##############")
# In this example, the solution by fitting is possible thanks to the prediction of water_ec without additional information.

sampleC4 = Soil(water =   [0.06,    0.08,   0.095,  0.11], 
                bulk_ec = [0.007,   0.0072, 0.0075, np.nan])

sample4EC = BulkEC(sampleC4)

print("BulkEC(sample4) :", sample4EC) # [0.007    0.0072   0.0075   0.007666]
print("sampleC4.df.water_ec :", sampleC4.df.water_ec)
print("sampleC4.Lw :", sampleC4.Lw)

#sampleC4.info.to_excel('sampleC4_info.xlsx')
#sampleC4.df.to_excel('sampleC4_df.xlsx')

sample4_Rhoades = Rhoades(sampleC4.df.water, 0.0264, 0.00628, 1.0, 0.38)

#fig = plt.figure()
#ax2e = fig.add_subplot(1, 1, 1)

#aa = 0.5
#ss = 50
#ax2e.set_title('Example 4')
#ax2e.scatter(sampleC4.df.water, sampleC4.bulk_ec, marker='D', color='black', s=ss)
#ax2e.plot(sampleC4.df.water, sample4EC, 'ro', alpha=aa, markersize=8)
#ax2e.plot(sampleC4.df.water, sample4_Rhoades, 'gD', alpha=aa, markersize=8)

#plt.show()


print("################## Example5 ###################")
        #                                       0       1       2      3       4     5
sampleC5 = Soil(water =                     [   0.06,   0.08,   0.095, 0.128             ], 
                            bulk_ec =       [   0.01,   0.014,  0.016, 0.02,   0.03, 0.04], 
                            temperature=25.+273.15)

print("BulkEC(sampleC5) :", BulkEC(sampleC5)) #  [0.01    0.014   0.016  0.02    0.03  0.04 ]
#sampleC5.info.to_excel('sampleC5_info.xlsx')
#sampleC5.df.to_excel('sampleC5_df.xlsx')


print(" ###############  Example6 #################")
sampleC6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.11]), 
                            temperature=25.+273.15,
                            clay = 20,
                            water_ec = 0.2, 
                            bulk_density=1.7,
                            frequency_ec = 60)

print("BulkEC(sampleC6) :", BulkEC(sampleC6)) # [0.00463  0.006493 0.007996 0.009588]
#sampleC6.info.to_excel('sampleC6_info.xlsx')
#sampleC6.df.to_excel('sampleC6_df.xlsx')


print(" ###############  Example6b #################")
sampleC6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            temperature=25.+273.15,
                            clay = 20,
                            water_ec = 0.2, 
                            bulk_density=1.7,
                            solid_ec=0.001)

print("BulkEC(sampleC6b) :", BulkEC(sampleC6b)) # [0.005271 0.007135 0.008637 0.012259]
#sampleC6b.info.to_excel('sampleC6b_info.xlsx')
#sampleC6b.df.to_excel('sampleC6b_df.xlsx')


print("############### Example7 #################")
sampleC7 = Soil(water =      [0.06,   0.08,   0.095, 0.128], 
                            temperature=25.+273.15, texture = "Clay", 
                            instrument = 'HydraProbe', 
                            orgm = 0.4)

print("BulkEC(sampleC7) :", BulkEC(sampleC7)) # [nan nan nan nan]
print('sampleC7.df.frequency_ec', sampleC7.df.frequency_ec)
#sampleC7.info.to_excel('sampleC7_info.xlsx')
#sampleC7.df.to_excel('sampleC7_df.xlsx')


print("############# Example8 ##################")
sampleC8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                            bulk_density=1.6, 
                            sand = 20, 
                            silt = 60, 
                            water_ec = 0.05, 
                            CEC = 20., 
                            frequency_perm = 50e6)

print("BulkEC(sampleC8) :", BulkEC(sampleC8)) # [0.004501 0.006082 0.007294 0.010038]
print('soil.Lw', sampleC8.Lw)
#sampleC8.info.to_excel('sampleC8_info.xlsx')
#sampleC8.df.to_excel('sampleC8_df.xlsx')


print("############## Example9b #####################")       
                         #                    0     1    2     3       4       5       6       7
sampleC9b = Soil(water =                     [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07], 
                            bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                            bulk_density=1.7,
                            instrument= 'GPR',
                            water_ec = 0.1,
                            clay = 10,
                            frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9])

                                            #[ 6.   11.  9.    6.675   5.651   nan     8.      8.5]
print("BulkPerm(sampleC9b) :", BulkPerm(sampleC9b)) 
#sampleC9b.info.to_excel('sampleC9b_info.xlsx')
#sampleC9b.df.to_excel('sampleC9b_df.xlsx')


##############
########################## Testing changing EC frequency, fitting and non-fitting ###############################
##############


print("#################  Example11  ####################")
sampleC11 = Soil(bulk_perm = [np.nan, 7],
                            frequency_ec = [1000, 1000, 1000, 1000, 1000, np.nan, 1000, 1000, 1000, 1000,   1000, 10], 
                            water =        [0.1,  0.1,  0.1,  0.1,  0.1,  0.1,    0.1,  0.1,  0.1,  np.nan, 0,    0.1],
                            bulk_density=1.5,
                            water_ec = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)

# [0.008988 0.008988 0.008988 0.008988 0.008988 0.008388 0.008988 0.008988 0.008988      nan 0.       0.008388]
print("BulkEC(sampleC11) :",BulkEC(sampleC11)) 
#sampleC11.info.to_excel('sampleC11_info.xlsx')
#sampleC11.df.to_excel('sampleC11_df.xlsx')


print("#################  Example12  ####################")
sampleC12 = Soil(bulk_perm = [np.nan, 7],
                            #               0    1    2     3     4     5      6      7       8       9         10   11
                            frequency_ec = [10,  20,  30,   55,   70,   66,    0,     np.nan, 48,     100,      5       ], 
                            water = [       0.1, 0.2, 0.11, 0.32, 0.61, 0.41,  0.01,  0.151,  0.21,   np.nan,   0,   0.1],
                            bulk_density=1.5,
                            water_ec = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
# [0.008388 0.017777 0.009282 0.030363 0.066724 0.04000794 0.013051 0.018771      nan 0.       0.008388]
print("BulkEC(sampleC12) :",BulkEC(sampleC12))
#sampleC12.info.to_excel('sampleC12_info.xlsx')
#sampleC12.df.to_excel('sampleC12_df.xlsx')


print("#################  Example13  ####################")
sampleC13 = Soil(bulk_perm = [np.nan, 7],
                            frequency_ec = [10,  10,  10,  10,  10,  np.nan, 10,  10,  10,  10,     10, 10], 
                            water =        [0.1, 0.1, 0.1, 0.1, 0.1, 0.1,    0.1, 0.1, 0.1, np.nan, 0,  0.1],
                            bulk_density=1.5,
                            water_ec = 0.05,
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            bulk_perm_inf = 5)
# [0.008388 0.008388 0.008388 0.008388 0.008388 0.00008388 0.008388 0.008388      nan 0.       0.008388]
print("BulkEC(sampleC13) :",BulkEC(sampleC13)) 
#sampleC13.info.to_excel('sampleC13_info.xlsx')
#sampleC13.df.to_excel('sampleC13_df.xlsx')


print("############ Example14 ##################")
sampleC14 = Soil(water =               [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                            bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                            frequency_ec = 10,
                            sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

# [0.0072   0.009    0.01          nan 0.0079   0.014    0.0079   0.010846 0.010467 0.006706 0.       0.347533]
print("BulkEC(sampleC14) :", BulkEC(sampleC14))
print("sampleC14.Lw :", sampleC14.Lw)
#sampleC14.info.to_excel('sampleC14_info.xlsx')
#sampleC14.df.to_excel('sampleC14_df.xlsx')


print("############ Example14b ##################")
sampleC14b = Soil(water =             [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                            bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                            frequency_ec = [ 1000],
                            sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
# [0.00774  0.009629 0.010676      nan 0.084163 0.0.084163 0.14033 0.133286 0.00722  0.       0.35574 ]
print("BulkEC(sampleC14b) :", BulkEC(sampleC14b))
print("sampleC14b.Lw :", sampleC14b.Lw)
#sampleC14b.info.to_excel('sampleC14b_info.xlsx')
#sampleC14b.df.to_excel('sampleC14b_df.xlsx')


print("############ Example14c ##################")
sampleC14c = Soil(water =                  [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2,    0.19, 0.01, 0, 0.5], 
                            bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014  ], 
                            frequency_ec = [ 10,     20,    30,   55,     70,     66,     0,    np.nan, 48,   99,   5 ],
                            sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

print("BulkEC(sample14c) :", BulkEC(sampleC14c))
print("sample14c.Lw :", sampleC14c.Lw)
# [0.0072   0.009    0.01          nan 0.0079   0.014    0.0079   0.010846 0.010467 0.006706 0.       0.347533]
#sampleC14c.info.to_excel('sampleC14c_info.xlsx')
#sampleC14c.df.to_excel('sampleC14c_df.xlsx')


print("############ Example14d ##################")
sampleC14d = Soil(water =                  [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0,   0.5], 
                            bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014                                  ], 
                            frequency_ec = [ 1000,   1500,  120,  150,    200,    500,    101,  100, 800,  1200, 1e6, 1e8],
                            sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

print("BulkEC(sampleC14d) :", BulkEC(sampleC14d))
print("sampleC14d.Lw :", sampleC14d.Lw)
#            [0.00774 0.009675 0.010498 nan  0.083438  0.014778 0.083156 0.138759 0.133112 0.007236 0.       0.446417]

#sampleC14d.info.to_excel('sampleC14d_info.xlsx')
#sampleC14d.df.to_excel('sampleC14d_df.xlsx')


print("############ Example14e ##################")
sampleC14e = Soil(water =              [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                            bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                            frequency_ec = 1000,
                            sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

print("BulkEC(sampleC14e) :", BulkEC(sampleC14e))
print("sampleC14e.Lw :", sampleC14e.Lw)
# [0.0072, 0.009, 0.01, nan, 0.081128, 0.014, 0.081128, 0.136013,            0.129118,           0.006706, 0.0, 0.347533]
# [0.00774 0.009629 0.010676      nan 0.084163 0.014853 0.084163 0.14033 0.133286 0.00722  0.       0.355739]

#sampleC14e.info.to_excel('sampleC14e_info.xlsx')
#sampleC14e.df.to_excel('sampleC14e_df.xlsx')


print("#################  Graph example  ####################")
sample11g = Soil(water = np.array([0.1]), 
                            frequency_ec=np.logspace(0, 9, 100), 
                            temperature=25.+273.15, 
                            water_ec = 0.05,
                            bulk_density=1.5,
                            sand = 20, 
                            silt = 60)
print(sample11g.bulk_ec)
print(sample11g.df.bulk_ec)

sample12g = Soil(water = np.array([0.2]), 
                            frequency_ec=np.logspace(0, 9, 100), 
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            CEC = 20., 
                            water_ec = 0.05,
                            bulk_density=1.5,
                            bulk_perm_inf = 5)

sample13g = Soil(water = np.array([0.3]), 
                            frequency_ec=np.logspace(0, 9, 100), 
                            temperature=25.+273.15, 
                            sand = 20, 
                            silt = 60, 
                            water_ec = 0.05,
                            bulk_density=1.5)

sample14g = Soil(water = np.array([0.4]), 
                            frequency_ec=np.logspace(0, 9, 100), 
                            temperature=25., 
                            sand = 20, 
                            silt = 60, 
                            water_ec = 0.05,
                            bulk_density=1.5,
                            CEC = 20., 
                            bulk_perm_inf = 5)


pred_ls11g = BulkEC(sample11g)
pred_ls12g = BulkEC(sample12g)
pred_ls13g = BulkEC(sample13g)
pred_ls14g = BulkEC(sample14g)

#plt.semilogx(sample11g.frequency_ec, pred_ls11g, 'bo', markersize=2)
#plt.semilogx(sample12g.frequency_ec, pred_ls12g, 'ro', markersize=2)
#plt.semilogx(sample13g.frequency_ec, pred_ls13g, 'go', markersize=2)
#plt.semilogx(sample14g.frequency_ec, pred_ls14g, 'yo', markersize=4)

#locs, labels = plt.xticks()  # Get the current locations and labels.
#plt.xticks(np.logspace(0, 8, 90), rotation='vertical')  # Set label locations.

#plt.show()