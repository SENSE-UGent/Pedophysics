import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
sys.path.insert(0, 'C:\\Users\\mendo\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')


from pedophysics.predict import BulkEC, BulkPerm, LongmireSmithP
from pedophysics.simulate import Soil
from pedophysics.utils.similar_arrays import arrays_are_similar

from pedophysics.pedophysical_models.bulk_ec import Rhoades

################################################################################################################
################################################## PREDICT BULK EC #############################################
################################################################################################################

########################## Testing EC DC frequency, fitting and non-fitting #########################


def test_sample_C0():
        sample_C0 = Soil(water = 0.1, 
                bulk_ec= [ 0.0072,    0.007,   0.0075,  0.008], 
                sand=20.0,
                silt = 40)

        assert (BulkEC(sample_C0) == np.array([0.0072, 0.007,  0.0075, 0.008])).all()


def test_sample_C0b():
        sample_C0b = Soil(water =              0.1, 
                bulk_ec = [ 0.0072,  0.007, 0.0075,  np.nan], 
                sand=20.0,
                silt = 10,
                bulk_density=1.5,
                water_ec = 0.05,
                instrument = 'GPR')
                                
        assert (BulkEC(sample_C0b) == np.array([0.0072, 0.007, 0.0075, 0.007 ])).all() 


def test_sample_C0c():
        sample_C0c = Soil(bulk_perm =                [np.nan, 7],
                        frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                        water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                        bulk_density=1.5,
                        water_ec = 0.05,
                        sand = 20, 
                        silt = 60, 
                        CEC = 20)

        expected_result = np.array([0.00866,  0.008765, 0.008815, 0.008867, 0.008924, 0.008988, np.nan, 0.008388, 0.009239, 0.009355, 0.009528, 0.009774])
        assert arrays_are_similar(BulkEC(sample_C0c), expected_result)


def test_sample_C0d():
        sample_C0d = Soil(bulk_perm =                [np.nan, 7],
                        frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                        water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                        bulk_density=1.5,
                        water_ec = 0.05,
                        texture = 'Silt loam',
                        instrument = 'EMI Dualem')

        expected_result = np.array([0.006533, 0.006611, 0.006654, 0.006691, 0.006739, 0.006795, np.nan, 0.006991, 0.007006, 0.007094, 0.007257, 0.007474])
        assert arrays_are_similar(BulkEC(sample_C0d), expected_result)


def test_sample_C1():
                                #                 0      1        2         3         4       5       6       7
        sample_C1 = Soil(water =                 [0.05,  0.1,     0.08,     0.11,     0.01,   np.nan, np.nan, 0.07    ], 
                        bulk_ec=         [0.006, 0.011,   0.009,    np.nan,   np.nan, np.nan, 0.008,  0.0085  ], 
                        water_ec = 0.1,
                        instrument = 'TDR')

        expected_result = np.array([0.006,  0.011, 0.009, 0.012123, np.nan, np.nan, 0.008, 0.0085])
        assert arrays_are_similar(BulkEC(sample_C1), expected_result)


def test_sample_C1b():
                                #                   0      1        2        3        4       5       6       7    
        sample_C1b = Soil(water =                  [0.05,  0.1,     0.08,    0.11,    0.01,   np.nan, np.nan, 0.07    ], 
                        bulk_ec =          [0.006, 0.011,   0.009,   np.nan,  np.nan, np.nan, 0.008,  0.0085  ], 
                        bulk_density=1.7,
                        water_ec = 0.1,
                        clay=2,
                        frequency_ec = 80)

        expected_result = np.array([0.006, 0.011, 0.009, 0.012277, 0.000123, np.nan, 0.008, 0.0085]) 
        assert arrays_are_similar(BulkEC(sample_C1b), expected_result)           


def test_sample_C1c():
        # Silt Loam sample Wunderlich et al., 2013 #### 0      1       2        3        4       5       6       7        8
        sample_C1c = Soil(water =                      [  0.06,  0.08,   0.1,     0.12,    0.14,   np.nan, 0.17,   0.19,    0.28    ], 
                        bulk_ec=         [  2e-3,  5.5e-2, np.nan,  np.nan,  np.nan, np.nan, 1.2e-2, np.nan,  np.nan  ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        water_ec = 0.1,
                        solid_perm = 5,
                        instrument = 'TDR')

        expected_result = np.array([0.002, 0.055, 0.009393, 0.012327, 0.01522, np.nan, 0.012, 0.022389, 0.011016])
        assert arrays_are_similar(BulkEC(sample_C1c), expected_result)   


def test_sample_C4():
# In this example, the solution by fitting is possible thanks to the prediction of water_ec without additional information.
        sample_C4 = Soil(water =   [0.06,    0.08,   0.095,  0.11], 
                        bulk_ec = [0.007,   0.0072, 0.0075, np.nan])

        expected_result = np.array([0.007, 0.0072, 0.0075, 0.007666])
        assert arrays_are_similar(BulkEC(sample_C4), expected_result)


def test_sample_C5():
                #                           0       1       2      3       4     5
        sample_C5 = Soil(water =        [   0.06,   0.08,   0.095, 0.128             ], 
                        bulk_ec =       [   0.01,   0.014,  0.016, 0.02,   0.03, 0.04], 
                        temperature=25.+273.15)

        expected_result = np.array([0.01, 0.014, 0.016, 0.02, 0.03, 0.04])
        assert arrays_are_similar(BulkEC(sample_C5), expected_result)     


def test_sample_C6():
        sample_C6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.11]), 
                        temperature=25.+273.15,
                        clay = 20,
                        water_ec = 0.2, 
                        bulk_density=1.7,
                        frequency_ec = 60)

        expected_result = np.array([0.00463,  0.006493, 0.007996, 0.009588])
        assert arrays_are_similar(BulkEC(sample_C6), expected_result)


def test_sample_C6b():
        sample_C6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                        temperature=25.+273.15,
                        clay = 20,
                        water_ec = 0.2, 
                        bulk_density=1.7,
                        solid_ec=0.001)

        expected_result = np.array([0.005271, 0.007135, 0.008637, 0.012259])
        assert arrays_are_similar(BulkEC(sample_C6b), expected_result)


def test_sample_C7():
        sample_C7 = Soil(water =      [0.06,   0.08,   0.095, 0.128], 
                        temperature=25.+273.15, texture = "Clay", 
                        instrument = 'HydraProbe', 
                        orgm = 0.4)

        expected_result = np.array([np.nan, np.nan, np.nan, np.nan])
        assert arrays_are_similar(BulkEC(sample_C7), expected_result)


def test_sample_C8():
        sample_C8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                        bulk_density=1.6, 
                        sand = 20, 
                        silt = 60, 
                        water_ec = 0.05, 
                        CEC = 20., 
                        frequency_perm = 50e6)

        expected_result = np.array([0.004501, 0.006082, 0.007294, 0.010038])
        assert arrays_are_similar(BulkEC(sample_C8), expected_result)


def test_sample_C9b():
                                #                 0     1    2     3       4       5       6       7
        sample_C9b = Soil(water =                [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07], 
                        bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        bulk_density=1.7,
                        instrument= 'GPR',
                        water_ec = 0.1,
                        clay = 10,
                        frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9])

        expected_result = np.array([6, 11, 9, 6.675, 5.651, np.nan, 8, 8.5])
        assert arrays_are_similar(BulkEC(sample_C9b), expected_result)     

########################## Testing changing EC frequency, fitting and non-fitting #######################


def test_sample_C11():
        sample_C11 = Soil(bulk_perm = [np.nan, 7],
                        frequency_ec = [1000, 1000, 1000, 1000, 1000, np.nan, 1000, 1000, 1000, 1000,   1000, 10], 
                        water =        [0.1,  0.1,  0.1,  0.1,  0.1,  0.1,    0.1,  0.1,  0.1,  np.nan, 0,    0.1],
                        bulk_density=1.5,
                        water_ec = 0.05,
                        temperature=25.+273.15, 
                        sand = 20, 
                        silt = 60, 
                        CEC = 20., 
                        bulk_perm_inf = 5)

        expected_result = np.array([0.008988, 0.008988, 0.008988, 0.008988, 0.008988, 0.008388, 0.008988, 0.008988, 0.008988, np.nan, 0., 0.008388])
        assert arrays_are_similar(BulkEC(sample_C11), expected_result)


def test_sample_C12():
        sample_C12 = Soil(bulk_perm = [np.nan, 7],
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
        expected_result = np.array([0.008388, 0.017777, 0.009282, 0.030363, 0.066724, 0.04000794, 0.013051, 0.018771, np.nan, 0., 0.008388])
        assert arrays_are_similar(BulkEC(sample_C12), expected_result)


def test_sample_C13():
        sample_C13 = Soil(bulk_perm = [np.nan, 7],
                        frequency_ec = [10,  10,  10,  10,  10,  np.nan, 10,  10,  10,  10,     10, 10], 
                        water =        [0.1, 0.1, 0.1, 0.1, 0.1, 0.1,    0.1, 0.1, 0.1, np.nan, 0,  0.1],
                        bulk_density=1.5, water_ec = 0.05, temperature=25.+273.15, sand = 20, silt = 60, CEC = 20., bulk_perm_inf = 5)
        expected_result = np.array([0.008388, 0.008388, 0.008388, 0.008388, 0.008388, 0.00008388, 0.008388, 0.008388, np.nan, 0., 0.008388])
        assert arrays_are_similar(BulkEC(sample_C13), expected_result)


def test_sample_C14():
        sample_C14 = Soil(water =         [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                        bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                        frequency_ec = 10,
                        sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

        expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.0079, 0.014, 0.0079, 0.010846, 0.010467, 0.006706, 0., 0.347533])
        assert arrays_are_similar(BulkEC(sample_C14), expected_result)


def test_sample_C14b():
        sample_C14b = Soil(water =        [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                        bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                        frequency_ec = [ 1000],
                        sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
        
        expected_result = np.array([0.00774, 0.009629, 0.010676, np.nan, 0.084163, 0., 0.084163, 0.14033, 0.133286, 0.00722, 0., 0.35574])
        assert arrays_are_similar(BulkEC(sample_C14b), expected_result)


def test_sample_C14c():
        sample_C14c = Soil(water =             [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2,    0.19, 0.01, 0, 0.5], 
                        bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014  ], 
                        frequency_ec = [ 10,     20,    30,   55,     70,     66,     0,    np.nan, 48,   99,   5 ],
                        sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

        expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.0079, 0.014, 0.0079, 0.010846, 0.010467, 0.006706, 0., 0.347533])
        assert arrays_are_similar(BulkEC(sample_C14c), expected_result)


def test_sample_C14d():
        sample_C14d = Soil(water =             [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0,   0.5], 
                        bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014                                  ], 
                        frequency_ec = [ 1000,   1500,  120,  150,    200,    500,    101,  100, 800,  1200, 1e6, 1e8],
                        sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')

        expected_result = np.array([0.00774, 0.009675, 0.010498, np.nan, 0.083438, 0.014778, 0.083156, 0.138759, 0.133112, 0.007236, 0., 0.446417])
        assert arrays_are_similar(BulkEC(sample_C14d), expected_result)


def test_sample_C14e():
        sample_C14e = Soil(water =        [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                        bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                        frequency_ec = 1000,
                        sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
        expected_result = np.array([0.00774, 0.009629, 0.010676, np.nan, 0.084163, 0.014853, 0.084163, 0.14033, 0.133286, 0.00722, 0., 0.355739])
        assert arrays_are_similar(BulkEC(sample_C14e), expected_result)


################################################################################################################
################################################## PREDICT BULK PERM ###########################################
################################################################################################################

########################## Testing fixed frequency, fitting and non-fitting ###############################

def test_sample_P0():
                                #   0      1    2     3       4       5       6       7
        sample_P0 = Soil( water =   [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07], 
                        bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        frequency_perm = [7e8, 7e8],
                        bulk_density=1.7,
                        solid_perm = 5,
                        instrument = 'HydraProbe',
                        salinity = 0.1)

        expected_result = np.array([6., 11., 9., np.nan, np.nan, np.nan, 8., 8.5])
        assert arrays_are_similar(BulkPerm(sample_P0), expected_result)


def test_sample_P1():
                                #    0     1    2     3       4       5       6       7
        sample_P1 = Soil( water =   [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07], 
                        bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        bulk_density=1.7,
                        solid_perm = 5,
                        instrument = 'TDR')

        expected_result = np.array([6., 11., 9., np.nan, 4.666, np.nan, 8., 8.5])
        assert arrays_are_similar(BulkPerm(sample_P1), expected_result)


def test_sample_P1b():
                                #    0     1    2     3       4       5       6       7
        sample_P1b = Soil(water =   [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07]    , 
                        bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        bulk_density=1.7,
                        clay=2,
                        solid_perm = 5,
                        instrument= 'GPR')
        expected_result = np.array([6., 11., 9., 10.653, 4.666, np.nan, 8., 8.5])
        assert arrays_are_similar(BulkPerm(sample_P1b), expected_result)


def test_sample_P3():
        sample_P3 = Soil( water = 0.1, 
                        bulk_perm = [7.2,    7,   7.5,     8], 
                        sand=20.0)
        expected_result = np.array([7.2, 7., 7.5, 8.])
        assert arrays_are_similar(BulkPerm(sample_P3), expected_result)


def test_sample_P3b():
        sample_P3b = Soil(water = 0.1, 
                        bulk_perm =       [7.2,  7,   7.5,   np.nan], 
                        sand = 20.0,
                        silt = 10,
                        instrument = 'GPR')
        expected_result = np.array([7.2, 7., 7.5, 7.])
        assert arrays_are_similar(BulkPerm(sample_P3b), expected_result)


def test_sample_P4():
        sample_P4 = Soil( water =   [0.06,   0.08,   0.095, 0.128], 
                        bulk_perm = [7,      7.2,    7.5 ,  np.nan], 
                        temperature=25.)
        expected_result = np.array([7., 7.2, 7.5,  np.nan])
        assert arrays_are_similar(BulkPerm(sample_P4), expected_result)


def test_sample_P6():
        sample_P6 = Soil(water = [0.06,   0.08,   0.095, 0.128], 
                        temperature=25., 
                        frequency_perm = 60e6)
        expected_result = np.array([np.nan, np.nan, np.nan, np.nan])
        assert arrays_are_similar(BulkPerm(sample_P6), expected_result)


def test_sample_P6b():
        sample_P6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                        temperature=25., texture = "Clay", 
                        instrument = 'HydraProbe', CEC = 2, bulk_density = 1.3, orgm = 0.4)
        expected_result = np.array([4.115, 4.817, 5.377, 6.707])
        assert arrays_are_similar(BulkPerm(sample_P6b), expected_result)


def test_sample_P6c():
        sample_P6c = Soil(water =  [0.06,  0.08,   0.095,  0.128], 
                        bulk_density=1.6, 
                        sand = 20, 
                        silt = 60, 
                        CEC = 20., 
                        frequency_perm = 50e6)
        expected_result = np.array([8.661, 10.462, 11.791, 14.661])
        assert arrays_are_similar(BulkPerm(sample_P6c), expected_result)


def test_sample_Pv():
                                 # 0        1       2       3        4        5        6        7      
        sample_Pv = Soil( water = [0.03,    0.08,   0.15,   0.20,    0.22,    0.07,    0.12,    0.18     ] , 
                        bulk_density=1.4,
                        texture = 'Sand',
                        solid_perm = 5,
                        CEC = 1.6,
                        frequency_perm = np.array([50e6]))
        expected_result = np.array([3.697, 5.315, 8.092, 10.447, 11.476, 4.967, 6.828, 9.468])
        assert arrays_are_similar(BulkPerm(sample_Pv), expected_result)


def test_sample_P7():
                                #        0     1    2     3       4       5       6       7
        sample_P7 = Soil(water =        [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07] , 
                        bulk_perm =     [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ] , 
                        bulk_density=1.7,
                        instrument= 'GPR',
                        frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9]  )
        expected_result = np.array([6., 11., 9., np.nan, np.nan, np.nan, 8., 8.5])
        assert arrays_are_similar(BulkPerm(sample_P7), expected_result)


def test_sample_P7b():
                                #         0     1    2     3       4       5       6       7
        sample_P7b = Soil(water =        [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07] , 
                        bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ] , 
                        bulk_density=1.7,
                        instrument= 'GPR',
                        water_ec = 0.2,
                        clay = 10,
                        frequency_perm =  [2e9,  2e9, 2e9,  2e9,    2e9,    5e9]  )
        expected_result = np.array([6., 11., 9., 6.785, 5.657, np.nan, 8., 8.5])
        assert arrays_are_similar(BulkPerm(sample_P7b), expected_result)


def test_sample_P8():
        sample_P8 = Soil(bulk_ec =               [0.0128, 0.0128, 0.0128, 0.0128, 0.0128, 0.0128], 
                                bulk_perm =      [np.nan, 7] ,
                                frequency_perm = [1e6 ,   1e6,    50e6,   100e6,  500e6,  1e9], 
                                #                [47.582  7.      14.049  12.271  8.888   8.134]
                                temperature=25., 
                                sand = 20, 
                                silt = 60, 
                                CEC = 20., 
                                bulk_perm_inf = 5)
        expected_result = np.array([47.582, 7., 14.049, 12.271, 8.888, 8.134])
        assert arrays_are_similar(BulkPerm(sample_P8), expected_result)
        assert arrays_are_similar(LongmireSmithP(sample_P8.bulk_ec, sample_P8.bulk_perm_inf, sample_P8.frequency_perm), np.array([47.58214722, 47.58214722, 14.04895663, 12.27116452, 8.88788941, 8.13446703]))
        #sampleP8.info.to_excel('sampleP8_info.xlsx')
        #sampleP8.df.to_excel('sampleP8_df.xlsx')


