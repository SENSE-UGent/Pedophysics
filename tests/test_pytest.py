import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
sys.path.insert(0, 'C:\\Users\\mendo\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')


from pedophysics.predict import BulkEC, BulkPerm, ParticleDensity, Salinity, WaterEC, Water
from pedophysics.simulate import Soil
from pedophysics.utils.similar_arrays import arrays_are_similar

from pedophysics.pedophysical_models.bulk_ec import Rhoades
from pedophysics.pedophysical_models.bulk_perm import LongmireSmithP

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
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      assert (BulkEC(sample_C0b) == np.array([0.0072, 0.007, 0.0075, 0.007])).all()


def test_sample_C0c():
      sample_C0c = Soil(bulk_perm =                [np.nan, 7],
            frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
            water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
            bulk_density=1.5, water_ec = 0.05, sand = 20, silt = 60, CEC = 20)
      expected_result = np.array([0.00866, 0.008765, 0.008815, 0.008867, 0.008924, 0.008988, np.nan, 0.008388, 0.009239, 0.009355, 0.009528, 0.009774])
      assert arrays_are_similar(BulkEC(sample_C0c), expected_result)


def test_sample_C0d():
        sample_C0d = Soil(bulk_perm =                [np.nan, 7],
                        frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                        water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                        bulk_density=1.5, water_ec = 0.05, texture = 'Silt loam', instrument = 'EMI Dualem')
        expected_result = np.array([0.006533, 0.006611, 0.006654, 0.006691, 0.006739, 0.006795, np.nan, 0.006991, 0.007006, 0.007094, 0.007257, 0.007474])
        assert arrays_are_similar(BulkEC(sample_C0d), expected_result)


def test_sample_C1():
                              #                 0      1        2         3         4       5       6       7
      sample_C1 = Soil(water =                 [0.05,  0.1,     0.08,     0.11,     0.01,   np.nan, np.nan, 0.07    ], 
                              bulk_ec=         [0.006, 0.011,   0.009,    np.nan,   np.nan, np.nan, 0.008,  0.0085  ], 
                              water_ec = 0.1,
                              instrument = 'TDR')
      expected_result = np.array([0.006,  0.011, 0.009, 0.012123, np.nan, np.nan, 0.008, 0.0085])
      print(BulkEC(sample_C1))
      assert arrays_are_similar(BulkEC(sample_C1), expected_result)


def test_sample_C1b():
                              #          0      1        2        3        4       5       6       7    
      sample_C1b = Soil(water =         [0.05,  0.1,     0.08,    0.11,    0.01,   np.nan, np.nan, 0.07    ], 
                              bulk_ec = [0.006, 0.011,   0.009,   np.nan,  np.nan, np.nan, 0.008,  0.0085  ], 
                              bulk_density=1.7, water_ec = 0.1, clay=2, frequency_ec = 80)
      expected_result = np.array([0.006, 0.011, 0.009, 0.012147, 0.000144, np.nan, 0.008, 0.0085]) 
      print(BulkEC(sample_C1b))
      assert arrays_are_similar(BulkEC(sample_C1b), expected_result)   


def test_sample_C1c():
      # Silt Loam sample Wunderlich et al., 2013 ####     0         1        2        3        4       5       6        7        8
      sample_C1c = Soil(water =        [  0.06,     0.12,    0.1,     0.12,    0.14,   np.nan, 0.23,    0.185,   0.28    ], 
                  bulk_ec=         [  1/(2e3),  1/(2e2), np.nan,  np.nan,  np.nan, np.nan, 1/(6e1), np.nan,  np.nan  ], 
                  bulk_density=1.7, texture = 'Sand', water_ec = 0.1, solid_perm = 5, instrument = 'TDR')
      bulkec_C1c = BulkEC(sample_C1c)
      print(bulkec_C1c)
      expected_result = np.array([0.0005, 0.005, 0.003713, 0.005565, 0.007463, np.nan, 0.01666667, 0.011825, 0.021218])
      assert arrays_are_similar(bulkec_C1c, expected_result) 


def test_sample_C4():
      # In this example, the solution by fitting is possible thanks to the prediction of water_ec without additional information.
      sample_C4 = Soil(water =   [0.06,    0.08,   0.095,  0.11], 
                  bulk_ec = [0.007,   0.0072, 0.0075, np.nan])
      BulkECsample_C4 = BulkEC(sample_C4)
      print(BulkECsample_C4 )
      expected_result = np.array([0.007, 0.0072, 0.0075, 0.007669])
      assert arrays_are_similar(BulkECsample_C4, expected_result)


def test_sample_C5():
                #                                   0       1       2      3       4     5
        sample_C5 = Soil(water =                [   0.06,   0.08,   0.095, 0.128             ], 
                                bulk_ec =       [   0.01,   0.014,  0.016, 0.02,   0.03, 0.04], 
                                temperature=25.+273.15)
        BulkECsample_C5 = BulkEC(sample_C5)
        expected_result = np.array([0.01, 0.014, 0.016, 0.02, 0.03, 0.04])
        assert arrays_are_similar(BulkECsample_C5, expected_result)


def test_sample_C6():
        sample_C6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.11]), 
                        temperature=25.+273.15, clay = 20, water_ec = 0.2, bulk_density=1.7, frequency_ec = 60)
        BulkECsample_C6 = BulkEC(sample_C6)
        expected_result = np.array([0.00489,  0.006819, 0.008372, 0.010014])
        assert arrays_are_similar(BulkECsample_C6, expected_result)


def test_sample_C6b():
      sample_C6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                  temperature=25.+273.15, clay = 20, water_ec = 0.2, bulk_density=1.7, solid_ec=0.001)
      BulkECsample_C6b = BulkEC(sample_C6b)
      print(BulkECsample_C6b)
      expected_result = np.array([0.005271, 0.007135, 0.008637, 0.012259])
      assert arrays_are_similar(BulkECsample_C6b, expected_result)


def test_sample_C7():
        sample_C7 = Soil(water =      [0.06,   0.08,   0.095, 0.128], 
                                temperature=25.+273.15, texture = "Clay", 
                                instrument = 'HydraProbe', 
                                orgm = 0.4)
        expected_result = np.array([np.nan, np.nan, np.nan, np.nan])
        assert arrays_are_similar(BulkEC(sample_C7), expected_result)


def test_sample_C8():
      sample_C8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                              bulk_density=1.6, sand = 20, silt = 60, water_ec = 0.05, CEC = 20., frequency_perm = 50e6)
      BulkECsample_C8 = BulkEC(sample_C8)
      print(BulkECsample_C8)
      expected_result = np.array([0.004501, 0.006082, 0.007294, 0.010038])
      assert arrays_are_similar(BulkECsample_C8, expected_result)


def test_sample_C9b():
                              #                 0     1    2     3       4       5       6       7
      sample_C9b = Soil(water =                [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07], 
                              bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                              bulk_density=1.7, instrument= 'GPR', water_ec = 0.1, clay = 10,
                              frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9])
      BulkECsample_C9b = BulkEC(sample_C9b)
      print(BulkECsample_C9b)
      expected_result = np.array([0.00188, 0.004261, 0.003249, 0.004797, 0.000336, np.nan, np.nan, 0.002772])
      assert arrays_are_similar(BulkECsample_C9b, expected_result)


def test_sample_C11():
      sample_C11 = Soil(bulk_perm = [np.nan, 7],
                              frequency_ec = [1000, 1000, 1000, 1000, 1000, np.nan, 1000, 1000, 1000, 1000,   1000, 10], 
                              water =        [0.1,  0.1,  0.1,  0.1,  0.1,  0.1,    0.1,  0.1,  0.1,  np.nan, 0,    0.1],
                              bulk_density=1.5, water_ec = 0.05, temperature=25.+273.15, sand = 20, silt = 60, CEC = 20., bulk_perm_inf = 5)
      BulkECsample_C11 = BulkEC(sample_C11)
      print(BulkECsample_C11)
      expected_result = np.array([0.008988, 0.008988, 0.008988, 0.008988, 0.008988, 0.008388, 0.008988, 0.008988, 0.008988, np.nan, 0., 0.00866])
      assert arrays_are_similar(BulkECsample_C11, expected_result)


def test_sample_C12():
      sample_C12 = Soil(bulk_perm = [np.nan, 7],
                              #               0    1    2     3     4     5      6      7       8       9         10   11
                              frequency_ec = [10,  20,  30,   55,   70,   66,    0,     np.nan, 48,     100,      5       ], 
                              water = [       0.1, 0.2, 0.11, 0.32, 0.61, 0.41,  0.01,  0.151,  0.21,   np.nan,   0,   0.1],
                              bulk_density=1.5, water_ec = 0.05, temperature=25.+273.15, sand = 20, silt = 60, CEC = 20., bulk_perm_inf = 5)
      BulkECsample_C12 = BulkEC(sample_C12)
      expected_result = np.array([0.00866, 0.018296, 0.009658, 0.031303, 0.068356, 0.041931, 0.000794, 0.013051, 0.019433, np.nan, 0., 0.008388])
      assert arrays_are_similar(BulkEC(sample_C12), expected_result)


def test_sample_C13():
      sample_C13 = Soil(bulk_perm = [np.nan, 7],
                              frequency_ec = [10,  10,  10,  10,  10,  np.nan, 10,  10,  10,  10,     10, 10], 
                              water =        [0.1, 0.1, 0.1, 0.1, 0.1, 0.1,    0.1, 0.1, 0.1, np.nan, 0,  0.1],
                              bulk_density=1.5,
                              water_ec = 0.05,
                              temperature=25.+273.15, 
                              sand = 20, 
                              silt = 60, 
                              CEC = 20., 
                              bulk_perm_inf = 5)
      BulkECsample_C13 = BulkEC(sample_C13)
      print(BulkECsample_C13)
      expected_result = np.array([0.00866, 0.00866, 0.00866, 0.00866, 0.00866, 0.008388, 0.00866, 0.00866, 0.00866, np.nan, 0., 0.00866])
      assert arrays_are_similar(BulkECsample_C13, expected_result)


def test_sample_C14():
      sample_C14 = Soil(water =         [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                              bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                              frequency_ec = 10,
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14 = BulkEC(sample_C14)
      print(BulkECsample_C14)
      expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.007888, 0.014, 0.007888, 00.010787, 0.010411, 0.006945, 0., 0.348232])
      assert arrays_are_similar(BulkECsample_C14, expected_result)


def test_sample_C14b():
      sample_C14b = Soil(water =        [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.25], 
                              bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                              frequency_ec = [ 1000],
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14b = BulkEC(sample_C14b)
      print(BulkECsample_C14b)
      expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.007882, 0.014, 0.007882, 0.010787, 0.01041, 0.00722, 0., 0.17567])
      assert arrays_are_similar(BulkECsample_C14b, expected_result)


def test_sample_C14c():
      sample_C14c = Soil(water =             [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2,    0.19, 0.01, 0, 0.3], 
                              bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014  ], 
                              frequency_ec = [ 10,     20,    30,   55,     70,     66,     0,    np.nan, 48,   99,   5 ],
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14c = BulkEC(sample_C14c)
      print(BulkECsample_C14c)
      expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.007981, 0.014, 0.007607, 0.010374, 0.010439, 0.007073, 0., 0.20552])
      assert arrays_are_similar(BulkECsample_C14c, expected_result)


def test_sample_C14d():
      sample_C14d = Soil(water = [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0,   0.3], 
            bulk_ec =          [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014                                  ], 
            frequency_ec =     [ 1000,   1500,  120,  150,    200,    500,    101,  100, 800,  1200, 1e6, 1e8],
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14d = BulkEC(sample_C14d)
      print(BulkECsample_C14d)
      expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.007799, 0.014, 0.007755, 0.010717, 0.010502, 0.007236, 0., 0.286932])
      assert arrays_are_similar(BulkECsample_C14d, expected_result)


def test_sample_C14e():
      sample_C14e = Soil(water =[ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.25 ], 
                  bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                  frequency_ec = 1000,
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14e = BulkEC(sample_C14e)
      print(BulkECsample_C14e)
      expected_result = np.array([0.0072, 0.009, 0.01, np.nan, 0.007882, 0.014, 0.007882, 0.010787, 0.01041, 0.00722, 0., 0.17567])
      assert arrays_are_similar(BulkECsample_C14e, expected_result)


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
                  temperature=25.+273.15, texture = "Clay", 
                  instrument = 'HydraProbe', CEC = 2, bulk_density = 1.3, orgm = 0.4)
      print(BulkPerm(sample_P6b))
      expected_result = np.array([4.09,  4.781, 5.331, 6.639])
      
      assert arrays_are_similar(BulkPerm(sample_P6b), expected_result)


def test_sample_P6c():
      sample_P6c = Soil(water =  [0.06,  0.08,   0.095,  0.128], 
                  bulk_density=1.6, 
                  sand = 20, 
                  silt = 60, 
                  CEC = 20., 
                  frequency_perm = 50e6)
      print(BulkPerm(sample_P6c))
      expected_result = np.array([ 8.531, 10.293, 11.593, 14.399])
      assert arrays_are_similar(BulkPerm(sample_P6c), expected_result)


def test_sample_Pv():
                                    # 0        1       2       3        4        5        6        7      
      sample_Pv = Soil( water = [0.03,    0.08,   0.15,   0.20,    0.22,    0.07,    0.12,    0.18     ] , 
                  bulk_density=1.4,
                  texture = 'Sand',
                  solid_perm = 5,
                  CEC = 1.6,
                  frequency_perm = np.array([50e6]))
      print(BulkPerm(sample_Pv))
      expected_result = np.array([ 3.687,  5.282,  8.014, 10.328, 11.339,  4.939,  6.771,  9.366])
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
      print(BulkPerm(sample_P7b))
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
      print(BulkPerm(sample_P8))
      expected_result = np.array([47.582, 7., 14.049, 12.271, 8.888, 8.134])
      assert arrays_are_similar(BulkPerm(sample_P8), expected_result)
      print(LongmireSmithP(sample_P8.bulk_ec, sample_P8.bulk_perm_inf, sample_P8.frequency_perm))
      assert arrays_are_similar(LongmireSmithP(sample_P8.bulk_ec, sample_P8.bulk_perm_inf, sample_P8.frequency_perm), np.array([47.58214722, 47.58214722, 14.04895663, 12.27116452, 8.88788941, 8.13446703]))
      #sampleP8.info.to_excel('sampleP8_info.xlsx')
      #sampleP8.df.to_excel('sampleP8_df.xlsx')

################################################################################################################
############################################### PREDICT PARTICLE DENSITY #######################################
################################################################################################################

def test_sample_PD1():
        sample_PD1 = Soil(water = np.append(np.random.rand(15)*40, np.full((4,1),np.nan)), 
                                bulk_perm = np.array([15, np.nan, 20, np.nan, 40]), 
                                instrument='TDR', 
                                particle_density=np.array([2, 2.2, 3, np.nan, 2.6]))

        expected_result = np.array([2., 2.2, 3., 2.65, 2.6, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65])
        assert arrays_are_similar(ParticleDensity(sample_PD1), expected_result)


def test_sample_PD2():
        sample_PD2 = Soil(water = np.append(np.random.rand(12)*40, np.full((4,1),np.nan)), 
                                sand = np.array([ np.nan, 30, 30,     20,     np.nan,  30,     np.nan, 20     ]), 
                                silt = np.array([ 10    , 30, np.nan, np.nan, np.nan,  np.nan, 20,     20     ]), 
                                clay = np.array([ np.nan, 30, 30,     np.nan, 20,      30,     20,     np.nan ]), 
                                orgm = np.array([ np.nan, 1,  np.nan, 1,      np.nan,  1,      0.5,    np.nan ]), 
                                particle_density=[2,      2,  2.2,    np.nan, np.nan,  np.nan, np.nan, np.nan] )

        expected_result = np.array([2., 2., 2.2, 2.65, 2.65, 2.69401972, 2.68255772, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65, 2.65])
        assert arrays_are_similar(ParticleDensity(sample_PD2), expected_result)

################################################################################################################
############################################### PREDICT SALINITY ###############################################
################################################################################################################


def test_sample_S1():
        sample_S1 = Soil( water_ec = np.array([100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3, )

        expected_result = np.array([0.00846, 0.01718, 0.00419, 0.02609, 0.00503, 0.00334, 0.0128, 0.02161, 0.01895, 0.02429, 0.02609, 0.04432])
        assert arrays_are_similar(Salinity(sample_S1), expected_result)


def test_sample_S2():
        sample_S2 = Soil( water_ec = np.array(  [100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3,
                        temperature = np.array([15,  10,  0,   40,   15, 15]) + 273 )

        expected_result = np.array([0.01089, 0.02589, 0.00937, 0.01966, 0.00647, 0.00429, 0.0128, 0.02161, 0.01895, 0.02429, 0.02609, 0.04432])
        assert arrays_are_similar(Salinity(sample_S2), expected_result)


def test_sample_Ss():
      sample_Ss = Soil( bulk_ec =  [0.02, 0.03, 0.04, 0.05, 0.06   ],
                  bulk_perm = [11.5, 14.8, 17,   20,   23    ],
                  clay=5,
                  bulk_density=1.48,
                  instrument = 'TDR')
      print(Salinity(sample_Ss))
      expected_result = np.array([0.02462, 0.02462, 0.02462, 0.02462, 0.02462])
      assert arrays_are_similar(Salinity(sample_Ss), expected_result)
      print(WaterEC(sample_Ss))
      assert arrays_are_similar(WaterEC(sample_Ss), np.array([0.283688, 0.283688, 0.283688, 0.283688, 0.283688]))


################################################################################################################
############################################### PREDICT WATEREC ################################################
################################################################################################################
# Testing Water EC from Bulk EC with examples from Brovelli & Cassiani 2011


def test_sample_ECW_DR_SCL():
        sample_ECW_DR_SCL = Soil( bulk_ec = [0, 1.6e-3, 4e-3, 9e-3, 1.5e-2, 2e-2],
                        water =  [0, 0.076,  0.15, 0.23, 0.3,    0.38])

        expected_result = np.array([0.068793, 0.068793, 0.068793, 0.068793, 0.068793, 0.068793])
        assert arrays_are_similar(WaterEC(sample_ECW_DR_SCL), expected_result)


def test_sample_ECW_DR_L():
        sample_ECW_DR_L = Soil( bulk_ec =  [0, 7*10**-3, 1.3*10**-2, 2*10**-2, 3*10**-2, 3.3*10**-2],
                        water = [0, 0.088,    0.18,       0.26,     0.35,     0.44      ])

        expected_result = np.array([0.09324, 0.09324, 0.09324, 0.09324, 0.09324, 0.09324])
        assert arrays_are_similar(WaterEC(sample_ECW_DR_L), expected_result)


def test_sample_ECW_DR_S():
        sample_ECW_DR_S = Soil( bulk_ec =  [0, 8*10**-4, 3*10**-3, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                        water = [0, 0.072,       0.144,       0.22,       0.29,       0.36])

        expected_result = np.array([0.063443, 0.063443, 0.063443, 0.063443, 0.063443, 0.063443])
        assert arrays_are_similar(WaterEC(sample_ECW_DR_S), expected_result)


def test_sample_ECW_DR_Sa():
        sample_ECW_DR_Sa = Soil( bulk_ec =  [0, 8*10**-4, np.nan, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                        water =  [0, 0.072,    0.144,  np.nan,     0.29,       0.36      ])

        expected_result = np.array([0.066926, 0.066926, 0.066926, 0.066926, 0.066926, 0.066926])
        assert arrays_are_similar(WaterEC(sample_ECW_DR_Sa), expected_result)


def test_sample_ECW_Odarslov_top():
        sample_ECW_Odarslov_top = Soil( bulk_ec = [0.02, 0.03, 0.04, 0.05, 0.06],
                        bulk_perm =           [11.5, 14.8,   17,   20,   23],
                        clay=5,
                        bulk_density=1.48,
                        instrument = 'TDR')

        expected_result = np.array([0.283688, 0.283688, 0.283688, 0.283688, 0.283688])
        assert arrays_are_similar(WaterEC(sample_ECW_Odarslov_top), expected_result)


def test_sample_ECW_Hil_ex():
        sample_ECW_Hil_ex = Soil( bulk_ec =        [0.025, 0.038, 0.065, 0.079, 0.1  ],
                        bulk_perm = [11.5,  15,    19,    22 ,   26   ],
                        clay=0,
                        bulk_density=1.8,
                        instrument = 'TDR')

        expected_result = np.array([0.427517, 0.427517, 0.427517, 0.427517, 0.427517])
        assert arrays_are_similar(WaterEC(sample_ECW_Hil_ex), expected_result)


def test_sample_ECW1():
        sample_ECW1 = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044])

        expected_result = np.array([0.09465033, 0.19791962, 0.04781566, 0.29905224, 0.05959418, 0.03598161,
        0.14085731, 0.24309688, 0.22055698, 0.27673502, 0.29905224, 0.49653988])
        assert arrays_are_similar(WaterEC(sample_ECW1), expected_result)


def test_sample_ECW2():
        sample_ECW2 = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044],
                        temperature = np.array([15,  10,  0,   40,   15, 15]) + 273)

        expected_result = np.array([0.07392045, 0.13297065, 0.02165263, 0.39356058, 0.04654703, 0.02810654,
        0.14085731, 0.24309688, 0.22055698, 0.27673502, 0.29905224, 0.49653988])
        assert arrays_are_similar(WaterEC(sample_ECW2), expected_result)   


################################################################################################################
####################################### PREDICT WATER FROM BULK PERM ###########################################
################################################################################################################


def test_sample_WP0():
        sample_WP0 = Soil( bulk_perm =        [10,   15, 20,   25,   7,  1,  12,  22,  5,  20,  30   ], 
                                bulk_density=1.7, texture = 'Sand', solid_perm = 5, instrument = 'GPR')

        expected_result = np.array([0.129, 0.21, 0.283, 0.35, 0.074, 0., 0.162, 0.31, 0.034, 0.283, 0.414])
        assert arrays_are_similar(Water(sample_WP0), expected_result)  


def test_sample_WP0b():
      sample_WP0b = Soil(water = [    0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                  bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ]), 
                  solid_perm = 5)

      expected_result = np.array([0.05, 0.11, 0.08, 0.11, np.nan, np.nan, np.nan, 0.07, np.nan, np.nan])
      assert arrays_are_similar(Water(sample_WP0b), expected_result)  


def test_sample_WP1():
      sample_WP1 = Soil(water =[0.05, 0.11, 0.08, np.nan, np.nan, 0.07      ], 
                  bulk_perm=   [6,    11,   9,    np.nan, 8,      8.5,    12], 
                  bulk_density=1.7,
                  instrument = 'TDR')

      expected_result = np.array([0.05, 0.11, 0.08, np.nan, 0.071, 0.07, 0.117])  
      assert arrays_are_similar(Water(sample_WP1), expected_result)  


def test_sample_WP1b():
      sample_WP1b = Soil(water =     [0.05, 0.11, 0.08, np.nan, np.nan, 0.07,   np.nan, 0.2,    0.02,   np.nan ], 
                  bulk_perm=         [6,    11,   9,    np.nan, 8,      8.5,    14,     np.nan, np.nan, 1      ], 
                  bulk_density=1.7,
                  texture = 'Sand',
                  solid_perm = 5,
                  instrument = 'TDR')

      expected_result = np.array([0.05, 0.11, 0.08, np.nan, 0.071, 0.07, 0.194, 0.2, 0.02, 0.])
      assert arrays_are_similar(Water(sample_WP1b), expected_result)  


def test_sample_WP1c():
      sample_WP1c = Soil(water =    [0.20, 0.31, 0.36, 0.38, 0.05                        ], 
                        bulk_perm=  [10,   15,   20,   25,   7,   1, 12, 22, 5, 20, 30   ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        solid_perm = 5)

      expected_result = np.array([0.2, 0.31, 0.36, 0.38, 0.05, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
      assert arrays_are_similar(Water(sample_WP1c), expected_result)  


def test_sample_WP3():
      sample_WP3 = Soil(water =     [0.20, 0.30, 0.35                                             ], 
                        bulk_perm = [10,   15,   20,   8.5, 8,    1,   12,   22,   5,    20,   30   ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        solid_perm = 5,
                        instrument = 'GPR')

      expected_result = np.array([0.2, 0.3, 0.35, 0.173, 0.164, 0., 0.234, 0.392, 0.104, 0.362, 0.414])
      assert arrays_are_similar(Water(sample_WP3), expected_result)  


def test_sample_WP4():
      sample_WP4 = Soil(water =            [0.20, 0.30, 0.35                                             ], 
                        bulk_perm=         [10,   15,   20,  8.5,   8,   1,   12,   22,   5,    20,   30 ], 
                        bulk_density=1.7,
                        clay = 40,
                        solid_perm = 5,
                        instrument = 'TDR')

      expected_result = np.array([0.2, 0.3, 0.35, 0.173, 0.164, 0., 0.234, 0.392, 0.104, 0.362, 0.459])
      assert arrays_are_similar(Water(sample_WP4), expected_result)  


def test_sample_WP5():
      sample_WP5 = Soil(water =     [0.20, 0.30, 0.35                              ], 
                        bulk_perm = [10,   15,   20,  8.5, 8, 1, 12, 22, 5, 20, 30 ], 
                        bulk_density=1.7, alpha = 0.3, frequency_perm = [150e6])

      expected_result = np.array([0.2, 0.3, 0.35, 0.173, 0.164, 0., 0.234, 0.392, 0.104, 0.362, 0.534])
      assert arrays_are_similar(Water(sample_WP5), expected_result)  


def test_sample_WP7b():
      sample_WP7b = Soil(water =     [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                  bulk_perm =        [6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5   ], 
                  frequency_perm =   [50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6 ],
                  bulk_density=1.7,
                  texture = 'Sand',
                  solid_perm = 5)

      expected_result = np.array([0.05,  0.11,  0.08,  0.11,  0.05,    np.nan, 0.056, 0.07,  0.061, 0.11 ]) # no water ec
      assert arrays_are_similar(Water(sample_WP7b), expected_result)  


#def test_sample_WP7c():
sample_WP7c = Soil(water =     [0.05,   0.11,   0.08,   0.11,   np.nan, np.nan, np.nan, 0.07,   np.nan, np.nan], 
            bulk_perm =        [6,      11,     9,      np.nan, 1,      np.nan, 8,      8.5,    8.5,    8.5   ],
            bulk_ec =          [np.nan, np.nan, np.nan, 0.002,  np.nan, np.nan, 0.003, np.nan, np.nan, np.nan],
            frequency_perm =   [50e6,   50e6,   50e6,   200e6,  200e6,  200e6,  50e6,   50e6,   50e6,   200e6 ],
            frequency_ec =     [50e2,   50e2,   50e2,   200e1,  200e3,  200e4,  50e3,   50e3,   50e3,   20    ],
            bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05)

expected_result = np.array([0.05, 0.11, 0.08, 0.11, 0., np.nan, 0.111, 0.07, 0., 0.117])
water_sample_WP7c = Water(sample_WP7c)
sample_WP7c.info.water
print(water_sample_WP7c)
#assert arrays_are_similar(water_sample_WP7c, expected_result)  


def test_sample_WP8():
      sample_WP8 = Soil( bulk_perm = [10,   15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
            bulk_density=1.7,
            texture = 'Sand',
            solid_perm = 5,
            frequency_perm =         [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6])

      expected_result = np.array([0.006, 0.012, 0.025, 0.046, 0.006, 0.005, 0.03, 0.107, np.nan, 0.25, 0.454])
      assert arrays_are_similar(Water(sample_WP8), expected_result)  


def test_sample_WP8b():
        sample_WP8b = Soil( bulk_perm =    [10,   15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
                        frequency_perm = [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6],
                        bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1)

        expected_result = np.array([0.019, 0.034, 0.081, 0.169, 0.018, 0.018, 0.105, 0.437, np.nan, 0.65, 0.65])
        assert arrays_are_similar(Water(sample_WP8b), expected_result)  


def test_sample_WP8c():
      sample_WP8c = Soil( bulk_perm =    [10,    15,    20,    25,    7,     1,    12,    np.nan, 5,      20,    22 ], 
                        frequency_perm = [1e6,   2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,   np.nan, 100e6, 200e6],
                        bulk_density=1.7, texture = 'Clay', solid_perm = 5, water_ec = 0.1)

      expected_result = np.array([0., 0.001, 0.003, 0.011, 0., 0., 0.005, np.nan, np.nan, 0.272, 0.65])
      assert arrays_are_similar(Water(sample_WP8c), expected_result) 


def test_sample_WP9():
      sample_WP9 = Soil( bulk_perm = [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
            bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_perm = 20e6)

      expected_result = np.array([0.044, 0.196, 0.376, 0.623, 0.019, 0.018, 0.085, 0.376, 0.018, 0.376, 0.457])
      assert arrays_are_similar(Water(sample_WP9), expected_result) 


def test_sample_WP9b():
      sample_WP9b = Soil( bulk_perm = [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
            bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_perm = 1e6)

      expected_result = np.array([0.019, 0.025, 0.039, 0.073, 0.018, 0.018, 0.019, 0.039, 0.018, 0.039, 0.055])
      assert arrays_are_similar(Water(sample_WP9b), expected_result) 


def test_sample_WPv():
      sample_WPv = Soil( bulk_perm = [3,    8,       15,    20,    22,    7,    12,    18     ], 
                        bulk_density=1.4, texture = 'Sand', solid_perm = 5, CEC = 1.6, frequency_perm = 50e6)

      expected_result = np.array([0.005, 0.15, 0.286, 0.363, 0.391, 0.126, 0.233, 0.333])
      assert arrays_are_similar(Water(sample_WPv), expected_result) 

################################################################################################################
####################################### PREDICT WATER FROM BULK EC #############################################
################################################################################################################


def test_sample_WEC1():
      sample_WEC1 = Soil( bulk_ec = np.array([10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                bulk_density=1.7,
                texture = 'Sand',
                water_ec = 0.1)

      expected_result = np.array([0.264, 0.335, 0.394, 0.447, 0.213, 0.065, 0.293, 0.394, 0.173, 0.394, 0.416])
      assert arrays_are_similar(Water(sample_WEC1), expected_result) 


def test_sample_WEC1b():
      sample_WEC1b = Soil( bulk_ec=np.array(  [10,   15,   20,   25,     7,    1,     12,    22,    5,    20,    30   ])*1e-3, 
                bulk_density=1.7, texture = 'Sand',
                water_ec = np.array( [ 0.05, 0.06, 0.07, np.nan, 0.01, 0.1]))

      expected_result = np.array([0.35, 0.413, 0.459, np.nan, 0.442, 0.065, np.nan, np.nan, np.nan, np.nan, np.nan])
      assert arrays_are_similar(Water(sample_WEC1b), expected_result) 


def test_sample_WEC2():
      sample_WEC2 = Soil(water = np.array([0.20, 0.31, 0.36, 0.38, 0.05                                   ]), 
                        bulk_ec= np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20,   30   ])*1e-3, 
                        bulk_density=1.7, texture = 'Sand')

      expected_result = np.array([0.2, 0.31, 0.36, 0.38, 0.05, 0., 0.186, 0.396, 0., 0.357, 0.539])
      assert arrays_are_similar(Water(sample_WEC2), expected_result) 


def test_sample_WEC3():
      sample_WEC3 = Soil(water = np.array([0.20, 0.31, 0.36, 0.38, 0.05                                 ]), 
                        bulk_ec=np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20, 30   ])*1e-3, 
                        bulk_density=1.7, water_ec=0.5, texture = 'Sand')

      expected_result = np.array([0.2, 0.31, 0.36, 0.38, 0.05, 0., 0.214, 0.373, 0., 0.349, 0.446])
      assert arrays_are_similar(Water(sample_WEC3), expected_result) 


def test_sample_WEC4():
      sample_WEC4 = Soil(water =    [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07], 
                  bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5 ])*1e-3, 
                  water_ec = 0.05)

      expected_result = np.array([0.05, 0.11, 0.08, 0.11, np.nan, np.nan, 0.067, 0.07])
      assert arrays_are_similar(Water(sample_WEC4), expected_result) 


def test_sample_WEC4b():
      sample_WEC4b = Soil(water = np.array([     0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                              bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', water_ec = 0.05)

      expected_result = np.array([0.05, 0.11, 0.08, 0.11, 0.079, np.nan, 0.067, 0.07, 0.073, 0.073])
      assert arrays_are_similar(Water(sample_WEC4b), expected_result) 


def test_sample_WEC5():
      sample_WEC5 = Soil( bulk_ec=np.array([          10,   15,   20,   25,   7,    1,    12,    20,  5,    20,    22 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', water_ec = 0.01, frequency_ec = 500)

      expected_result = np.array([0.566, 0.65, 0.65, 0.65, 0.421, 0.15, 0.639, 0.65, 0.316, 0.65, 0.65 ])
      assert arrays_are_similar(Water(sample_WEC5), expected_result) 


def test_sample_WEC5b():
      sample_WEC5b = Soil( bulk_ec=np.array([           10,    0,    np.nan, np.nan,   7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_ec = np.array([2e6]))

      expected_result = np.array([0.216, 0.018, np.nan, np.nan, 0.168, 0.029, 0.242, 0.34, 0.13, 0.34, 0.361])
      assert arrays_are_similar(Water(sample_WEC5b), expected_result) 


def test_sample_WEC6():
      sample_WEC6 = Soil(water = np.array([0.20, 0.30, 0.35]), 
            bulk_ec=np.array( [10,   15,   20,   8.5,  8,    1,    12,   22,   5,    20,   30   ])*1e-3, 
            bulk_density=1.7, clay = 40, water_ec=0.4, frequency_ec=5e3)

      expected_result = np.array([0.2, 0.3, 0.35, 0.162, 0.148, 0.005, 0.241, 0.374, 0.033, 0.354, 0.131])
      assert arrays_are_similar(Water(sample_WEC6), expected_result) 


def test_sample_WEC6b():
      sample_WEC6b = Soil( bulk_ec=np.array([   10,   15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                                    water =         [0.1,  0.12], bulk_density=1.7, texture = 'Sand', water_ec = 0.1,
                        frequency_ec=np.array([1,    2,     2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

      print('sample_WEC6b.Lw', sample_WEC6b.Lw)
      expected_result = np.array([0.1, 0.12, 0.394, 0.447, 0.213, 0.063, 0.287, 0.388, 0.168, 0.386, 0.406])
      assert arrays_are_similar(Water(sample_WEC6b), expected_result)    


def test_sample_WEC6c():
      sample_WEC6c = Soil( bulk_ec=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
            frequency_ec = np.array([1,    2,      2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]),
            bulk_density=1.7, texture = 'Clay', water_ec = 0.1)

      expected_result = np.array([0.023, 0.034, 0.045, 0.056, 0.016, 0.002, 0.026, 0.044, 0.011, 0.043, 0.047])
      assert arrays_are_similar(Water(sample_WEC6c), expected_result)  


def test_sample_WEC7():
      sample_WEC7 = Soil(water =           [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                        bulk_ec =np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                        bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05,
                        frequency_ec =    [50e1, 50e2, 50e2, 200e2,  200e2,  200e2,  50e2,   50e1, 50e1,   200e2] )
      
      expected_result = np.array([0.05, 0.11, 0.08, 0.11, 0.071, np.nan, 0.028, 0.07, 0.028, 0.028])
      assert arrays_are_similar(Water(sample_WEC7), expected_result)  


def test_sample_WEC7b():
      sample_WEC7b = Soil(water = np.array(    [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                              bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05,
                        frequency_ec =         [50,   5,   50,    200e2,  200e2,  200e2,  50e2,   50e1, 50,     20])

      expected_result = np.array([0.05, 0.11, 0.08, 0.11, 0.071, np.nan, 0.063, 0.07, 0.072, 0.073])
      assert arrays_are_similar(Water(sample_WEC7b), expected_result)  


def test_sample_WECv():
      sample_WECv = Soil( bulk_ec=np.array([3,    8,    15,   20,   22,    7,     12,    18,   10,   2     ])*1e-3, 
                  frequency_ec = np.array( [50,   500,  5000, 2000, 50000, 55000, 25000, 100,  10,   50]),
                  bulk_density=1.4, texture = 'Sand', CEC = 1.6, water_ec = 0.1)

      expected_result = np.array([0.111, 0.209, 0.306, 0.366, 0.378, 0.18, 0.262, 0.348, 0.245, 0.083])
      assert arrays_are_similar(Water(sample_WECv), expected_result) 



