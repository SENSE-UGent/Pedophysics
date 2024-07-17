import sys
import os
import pandas as pd
import ast

# Get notebook and parent dir
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)

# Set path to pedophysics module 
pedophysics_code_path = os.path.join(parent_dir)
sys.path.insert(0, pedophysics_code_path)

import numpy as np

from pedophysics.predict import BulkEC, BulkPerm, ParticleDensity, Salinity, WaterEC, Water
from pedophysics.simulate import Soil
from pedophysics.utils.similar_arrays import arrays_are_similar

from pedophysics.pedophysical_models.bulk_ec import Rhoades
from pedophysics.pedophysical_models.bulk_perm import LongmireSmithP

############################################# LOAD TEST DATA ############################################

def load_expected_results(csv_file):
    parent_dir_path = 'tests/'  # This navigates one directory up from the current folder
    full_path = parent_dir_path + csv_file
    df = pd.read_csv(full_path)
    safe_env = {'nan': np.nan}  # Define a safe environment with only np.nan
    return {
        col: np.array(eval(df[col][0], {"__builtins__": None}, safe_env))
        for col in df.columns
    }

# Load expected results
expected_results = load_expected_results('test_data.csv')

################################################################################################################
################################################## PREDICT BULK EC #############################################
################################################################################################################

########################## Testing EC DC frequency, fitting and non-fitting #########################


def test_sample_C0():
      sample_C0 = Soil(water = 0.1, 
                  bulk_ec= [ 0.0072,    0.007,   0.0075,  0.008], 
                  sand=20.0,
                  silt = 40)
      assert arrays_are_similar(BulkEC(sample_C0), expected_results['test_sample_C0'])


def test_sample_C0b():
      sample_C0b = Soil(water =              0.1, 
                  bulk_ec = [ 0.0072,  0.007, 0.0075,  np.nan], 
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      assert arrays_are_similar(BulkEC(sample_C0b), expected_results['test_sample_C0b'])


def test_sample_C0c():
      sample_C0c = Soil(bulk_perm =                [np.nan, 7],
            frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
            water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
            bulk_density=1.5, water_ec = 0.05, sand = 20, silt = 60, CEC = 20)
      
      BulkEC_C0c = BulkEC(sample_C0c)
      print('BulkEC_C0c', BulkEC_C0c)
      assert arrays_are_similar(BulkEC_C0c, expected_results['test_sample_C0c'])


def test_sample_C0d():
        sample_C0d = Soil(bulk_perm =                [np.nan, 7],
                        frequency_ec = [10 ,     50,      100,     200,     500,     1000,     2000,   np.nan,  10000,   20000,   50000,   1e5], 
                        water =        [0.1,     0.1,     0.1,     0.1,     0.1,     0.1,      np.nan, 0.1,     0.1,     0.1,     0.1,     0.1],
                        bulk_density=1.5, water_ec = 0.05, texture = 'Silt loam', instrument = 'EMI Dualem')
        assert arrays_are_similar(BulkEC(sample_C0d), expected_results['test_sample_C0d'])


def test_sample_C1():
                              #                 0      1        2         3         4       5       6       7
      sample_C1 = Soil(water =                 [0.05,  0.1,     0.08,     0.11,     0.01,   np.nan, np.nan, 0.07    ], 
                              bulk_ec=         [0.006, 0.011,   0.009,    np.nan,   np.nan, np.nan, 0.008,  0.0085  ], 
                              water_ec = 0.1,
                              instrument = 'TDR')
      print(BulkEC(sample_C1))
      assert arrays_are_similar(BulkEC(sample_C1), expected_results['test_sample_C1'])


def test_sample_C1b():
                              #          0      1        2        3        4       5       6       7    
      sample_C1b = Soil(water =         [0.05,  0.1,     0.08,    0.11,    0.01,   np.nan, np.nan, 0.07    ], 
                              bulk_ec = [0.006, 0.011,   0.009,   np.nan,  np.nan, np.nan, 0.008,  0.0085  ], 
                              bulk_density=1.7, water_ec = 0.1, clay=2, frequency_ec = 80)
      print(BulkEC(sample_C1b))
      assert arrays_are_similar(BulkEC(sample_C1b), expected_results['test_sample_C1b'])   


def test_sample_C1c():
      # Silt Loam sample Wunderlich et al., 2013 ####     0         1        2        3        4       5       6        7        8
      sample_C1c = Soil(water =        [  0.06,     0.12,    0.1,     0.12,    0.14,   np.nan, 0.23,    0.185,   0.28    ], 
                  bulk_ec=         [  1/(2e3),  1/(2e2), np.nan,  np.nan,  np.nan, np.nan, 1/(6e1), np.nan,  np.nan  ], 
                  bulk_density=1.7, texture = 'Sand', water_ec = 0.1, solid_perm = 5, instrument = 'TDR')
      bulkec_C1c = BulkEC(sample_C1c)
      print(bulkec_C1c)

      assert arrays_are_similar(bulkec_C1c, expected_results['test_sample_C1c']) 


def test_sample_C4():
      # In this example, the solution by fitting is possible thanks to the prediction of water_ec without additional information.
      sample_C4 = Soil(water =   [0.06,    0.08,   0.095,  0.11], 
                  bulk_ec = [0.007,   0.0072, 0.0075, np.nan])
      BulkECsample_C4 = BulkEC(sample_C4)
      print(BulkECsample_C4 )
      assert arrays_are_similar(BulkECsample_C4, expected_results['test_sample_C4'])


def test_sample_C5():
                #                                   0       1       2      3       4     5
        sample_C5 = Soil(water =                [   0.06,   0.08,   0.095, 0.128             ], 
                                bulk_ec =       [   0.01,   0.014,  0.016, 0.02,   0.03, 0.04], 
                                temperature=25.+273.15)
        BulkECsample_C5 = BulkEC(sample_C5)
        assert arrays_are_similar(BulkECsample_C5, expected_results['test_sample_C5'])


def test_sample_C6():
        sample_C6 = Soil(water = np.array([0.06,   0.08,   0.095, 0.11]), 
                        temperature=25.+273.15, clay = 20, water_ec = 0.2, bulk_density=1.7, frequency_ec = 60)
        BulkECsample_C6 = BulkEC(sample_C6)
        assert arrays_are_similar(BulkECsample_C6, expected_results['test_sample_C6'])


def test_sample_C6b():
      sample_C6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                  temperature=25.+273.15, clay = 20, water_ec = 0.2, bulk_density=1.7, solid_ec=0.001)
      BulkECsample_C6b = BulkEC(sample_C6b)
      print(BulkECsample_C6b)
      assert arrays_are_similar(BulkECsample_C6b, expected_results['test_sample_C6b'])


def test_sample_C7():
        sample_C7 = Soil(water =      [0.06,   0.08,   0.095, 0.128], 
                                temperature=25.+273.15, texture = "Clay", 
                                instrument = 'HydraProbe', 
                                orgm = 0.4)
        assert arrays_are_similar(BulkEC(sample_C7), expected_results['test_sample_C7'])


def test_sample_C8():
      sample_C8 = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                              bulk_density=1.6, sand = 20, silt = 60, water_ec = 0.05, CEC = 20., frequency_perm = 50e6)
      BulkECsample_C8 = BulkEC(sample_C8)
      print(BulkECsample_C8)
      assert arrays_are_similar(BulkECsample_C8, expected_results['test_sample_C8'])


def test_sample_C9b():
                              #                 0     1    2     3       4       5       6       7
      sample_C9b = Soil(water =                [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07], 
                              bulk_perm =      [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                              bulk_density=1.7, instrument= 'GPR', water_ec = 0.1, clay = 10,
                              frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9])
      BulkECsample_C9b = BulkEC(sample_C9b)
      print(BulkECsample_C9b)
      assert arrays_are_similar(BulkECsample_C9b, expected_results['test_sample_C9b'])


def test_sample_C11():
      sample_C11 = Soil(bulk_perm = [np.nan, 7],
                              frequency_ec = [1000, 1000, 1000, 1000, 1000, np.nan, 1000, 1000, 1000, 1000,   1000, 10], 
                              water =        [0.1,  0.1,  0.1,  0.1,  0.1,  0.1,    0.1,  0.1,  0.1,  np.nan, 0,    0.1],
                              bulk_density=1.5, water_ec = 0.05, temperature=25.+273.15, sand = 20, silt = 60, CEC = 20., bulk_perm_inf = 5)
      BulkECsample_C11 = BulkEC(sample_C11)
      print(BulkECsample_C11)
      assert arrays_are_similar(BulkECsample_C11, expected_results['test_sample_C11'])


def test_sample_C12():
      sample_C12 = Soil(bulk_perm = [np.nan, 7],
                              #               0    1    2     3     4     5      6      7       8       9         10   11
                              frequency_ec = [10,  20,  30,   55,   70,   66,    0,     np.nan, 48,     100,      5       ], 
                              water = [       0.1, 0.2, 0.11, 0.32, 0.61, 0.41,  0.01,  0.151,  0.21,   np.nan,   0,   0.1],
                              bulk_density=1.5, water_ec = 0.05, temperature=25.+273.15, sand = 20, silt = 60, CEC = 20., bulk_perm_inf = 5)
      BulkECsample_C12 = BulkEC(sample_C12)
      assert arrays_are_similar(BulkEC(sample_C12), expected_results['test_sample_C12'])


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
      assert arrays_are_similar(BulkECsample_C13, expected_results['test_sample_C13'])


def test_sample_C14():
      sample_C14 = Soil(water =         [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.5 ], 
                              bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                              frequency_ec = 10,
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14 = BulkEC(sample_C14)
      print(BulkECsample_C14)
      assert arrays_are_similar(BulkECsample_C14, expected_results['test_sample_C14'])


def test_sample_C14b():
      sample_C14b = Soil(water =        [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.25], 
                              bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                              frequency_ec = [ 1000],
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14b = BulkEC(sample_C14b)
      print(BulkECsample_C14b)
      assert arrays_are_similar(BulkECsample_C14b, expected_results['test_sample_C14b'])


def test_sample_C14c():
      sample_C14c = Soil(water =             [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2,    0.19, 0.01, 0, 0.3], 
                              bulk_ec =      [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014  ], 
                              frequency_ec = [ 10,     20,    30,   55,     70,     66,     0,    np.nan, 48,   99,   5 ],
                              sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14c = BulkEC(sample_C14c)
      print(BulkECsample_C14c)
      assert arrays_are_similar(BulkECsample_C14c, expected_results['test_sample_C14c'])


def test_sample_C14d():
      sample_C14d = Soil(water = [ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0,   0.3], 
            bulk_ec =          [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014                                  ], 
            frequency_ec =     [ 1000,   1500,  120,  150,    200,    500,    101,  100, 800,  1200, 1e6, 1e8],
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14d = BulkEC(sample_C14d)
      print(BulkECsample_C14d)
      assert arrays_are_similar(BulkECsample_C14d, expected_results['test_sample_C14d'])


def test_sample_C14e():
      sample_C14e = Soil(water =[ 0.1,    0.15,  0.18, np.nan, 0.12,   np.nan, 0.12, 0.2, 0.19, 0.01, 0, 0.25 ], 
                  bulk_ec = [ 0.0072, 0.009, 0.01, np.nan, np.nan, 0.014,  ], 
                  frequency_ec = 1000,
                  sand=20.0, silt = 10, bulk_density=1.5, water_ec = 0.05, instrument = 'GPR')
      BulkECsample_C14e = BulkEC(sample_C14e)
      print(BulkECsample_C14e)
      assert arrays_are_similar(BulkECsample_C14e, expected_results['test_sample_C14e'])


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

        assert arrays_are_similar(BulkPerm(sample_P0), expected_results['test_sample_P0'])


def test_sample_P1():
                                #    0     1    2     3       4       5       6       7
        sample_P1 = Soil( water =   [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07], 
                        bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        bulk_density=1.7,
                        solid_perm = 5,
                        instrument = 'TDR')

        assert arrays_are_similar(BulkPerm(sample_P1), expected_results['test_sample_P1'])


def test_sample_P1b():
                                #    0     1    2     3       4       5       6       7
        sample_P1b = Soil(water =   [0.05, 0.1, 0.08, 0.14,   0.04,   np.nan, np.nan, 0.07]    , 
                        bulk_perm = [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ], 
                        bulk_density=1.7,
                        clay=2,
                        solid_perm = 5,
                        instrument= 'GPR')
        
        BulkPerm_P1b = BulkPerm(sample_P1b)
        print('BulkPerm_P1b', BulkPerm_P1b)
        assert arrays_are_similar(BulkPerm_P1b, expected_results['test_sample_P1b'])


def test_sample_P3():
        sample_P3 = Soil( water = 0.1, 
                        bulk_perm = [7.2,    7,   7.5,     8], 
                        sand=20.0)
        assert arrays_are_similar(BulkPerm(sample_P3), expected_results['test_sample_P3'])


def test_sample_P3b():
        sample_P3b = Soil(water = 0.1, 
                        bulk_perm =       [7.2,  7,   7.5,   np.nan], 
                        sand = 20.0,
                        silt = 10,
                        instrument = 'GPR')
        assert arrays_are_similar(BulkPerm(sample_P3b), expected_results['test_sample_P3b'])


def test_sample_P4():
        sample_P4 = Soil( water =   [0.06,   0.08,   0.095, 0.128], 
                        bulk_perm = [7,      7.2,    7.5 ,  np.nan], 
                        temperature=25.)
        assert arrays_are_similar(BulkPerm(sample_P4), expected_results['test_sample_P4'])


def test_sample_P6():
        sample_P6 = Soil(water = [0.06,   0.08,   0.095, 0.128], 
                        temperature=25., 
                        frequency_perm = 60e6)
        assert arrays_are_similar(BulkPerm(sample_P6), expected_results['test_sample_P6'])


def test_sample_P6b():
      sample_P6b = Soil(water = np.array([0.06,   0.08,   0.095, 0.128]), 
                  temperature=25.+273.15, texture = "Clay", 
                  instrument = 'HydraProbe', CEC = 2, bulk_density = 1.3, orgm = 0.4)
      print(BulkPerm(sample_P6b))
      
      assert arrays_are_similar(BulkPerm(sample_P6b), expected_results['test_sample_P6b'])


def test_sample_P6c():
      sample_P6c = Soil(water =  [0.06,  0.08,   0.095,  0.128], 
                  bulk_density=1.6, 
                  sand = 20, 
                  silt = 60, 
                  CEC = 20., 
                  frequency_perm = 50e6)
      print(BulkPerm(sample_P6c))
      assert arrays_are_similar(BulkPerm(sample_P6c), expected_results['test_sample_P6c'])


def test_sample_Pv():
                                    # 0        1       2       3        4        5        6        7      
      sample_Pv = Soil( water = [0.03,    0.08,   0.15,   0.20,    0.22,    0.07,    0.12,    0.18     ] , 
                  bulk_density=1.4,
                  texture = 'Sand',
                  solid_perm = 5,
                  CEC = 1.6,
                  frequency_perm = np.array([50e6]))
      print(BulkPerm(sample_Pv))
      assert arrays_are_similar(BulkPerm(sample_Pv), expected_results['test_sample_Pv'])


def test_sample_P7():
                                #        0     1    2     3       4       5       6       7
        sample_P7 = Soil(water =        [0.05, 0.1, 0.08, 0.11,   0.01,   np.nan, np.nan, 0.07] , 
                        bulk_perm =     [6,    11,  9,    np.nan, np.nan, np.nan, 8,      8.5 ] , 
                        bulk_density=1.7,
                        instrument= 'GPR',
                        frequency_perm = [2e9,  2e9, 2e9,  2e9,    2e9,    5e9]  )
        assert arrays_are_similar(BulkPerm(sample_P7), expected_results['test_sample_P7'])


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
      assert arrays_are_similar(BulkPerm(sample_P7b), expected_results['test_sample_P7b'])


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
      assert arrays_are_similar(BulkPerm(sample_P8), expected_results['test_sample_P8'])
      #print(LongmireSmithP(sample_P8.bulk_ec, sample_P8.bulk_perm_inf, sample_P8.frequency_perm))
      #assert arrays_are_similar(LongmireSmithP(sample_P8.bulk_ec, sample_P8.bulk_perm_inf, sample_P8.frequency_perm), np.array([47.58214722, 47.58214722, 14.04895663, 12.27116452, 8.88788941, 8.13446703]))
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

        assert arrays_are_similar(ParticleDensity(sample_PD1), expected_results['test_sample_PD1'])


def test_sample_PD2():
        sample_PD2 = Soil(water = np.append(np.random.rand(12)*40, np.full((4,1),np.nan)), 
                                sand = np.array([ np.nan, 30, 30,     20,     np.nan,  30,     np.nan, 20     ]), 
                                silt = np.array([ 10    , 30, np.nan, np.nan, np.nan,  np.nan, 20,     20     ]), 
                                clay = np.array([ np.nan, 30, 30,     np.nan, 20,      30,     20,     np.nan ]), 
                                orgm = np.array([ np.nan, 1,  np.nan, 1,      np.nan,  1,      0.5,    np.nan ]), 
                                particle_density=[2,      2,  2.2,    np.nan, np.nan,  np.nan, np.nan, np.nan] )

        assert arrays_are_similar(ParticleDensity(sample_PD2), expected_results['test_sample_PD2'])

################################################################################################################
############################################### PREDICT SALINITY ###############################################
################################################################################################################


def test_sample_S1():
        sample_S1 = Soil( water_ec = np.array([100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3, )

        assert arrays_are_similar(Salinity(sample_S1), expected_results['test_sample_S1'])


def test_sample_S2():
        sample_S2 = Soil( water_ec = np.array(  [100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3,
                        temperature = np.array([15,  10,  0,   40,   15, 15]) + 273 )

        assert arrays_are_similar(Salinity(sample_S2), expected_results['test_sample_S2'])


def test_sample_Ss():
      sample_Ss = Soil( bulk_ec =  [0.02, 0.03, 0.04, 0.05, 0.06   ],
                  bulk_perm = [11.5, 14.8, 17,   20,   23    ],
                  clay=5,
                  bulk_density=1.48,
                  instrument = 'TDR')
      print(Salinity(sample_Ss))
      assert arrays_are_similar(Salinity(sample_Ss), expected_results['test_sample_Ss'])
      #print(WaterEC(sample_Ss))
      #assert arrays_are_similar(WaterEC(sample_Ss), np.array([0.283688, 0.283688, 0.283688, 0.283688, 0.283688]))


################################################################################################################
############################################### PREDICT WATEREC ################################################
################################################################################################################
# Testing Water EC from Bulk EC with examples from Brovelli & Cassiani 2011


def test_sample_ECW_DR_SCL():
        sample_ECW_DR_SCL = Soil( bulk_ec = [0, 1.6e-3, 4e-3, 9e-3, 1.5e-2, 2e-2],
                        water =  [0, 0.076,  0.15, 0.23, 0.3,    0.38])

        assert arrays_are_similar(WaterEC(sample_ECW_DR_SCL), expected_results['test_sample_ECW_DR_SCL'])


def test_sample_ECW_DR_L():
        sample_ECW_DR_L = Soil( bulk_ec =  [0, 7*10**-3, 1.3*10**-2, 2*10**-2, 3*10**-2, 3.3*10**-2],
                        water = [0, 0.088,    0.18,       0.26,     0.35,     0.44      ])

        assert arrays_are_similar(WaterEC(sample_ECW_DR_L), expected_results['test_sample_ECW_DR_L'])


def test_sample_ECW_DR_S():
        sample_ECW_DR_S = Soil( bulk_ec =  [0, 8*10**-4, 3*10**-3, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                        water = [0, 0.072,       0.144,       0.22,       0.29,       0.36])

        assert arrays_are_similar(WaterEC(sample_ECW_DR_S), expected_results['test_sample_ECW_DR_S'])


def test_sample_ECW_DR_Sa():
        sample_ECW_DR_Sa = Soil( bulk_ec =  [0, 8*10**-4, np.nan, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                        water =  [0, 0.072,    0.144,  np.nan,     0.29,       0.36      ])

        assert arrays_are_similar(WaterEC(sample_ECW_DR_Sa), expected_results['test_sample_ECW_DR_Sa'])


def test_sample_ECW_Odarslov_top():
        sample_ECW_Odarslov_top = Soil( bulk_ec = [0.02, 0.03, 0.04, 0.05, 0.06],
                        bulk_perm =           [11.5, 14.8,   17,   20,   23],
                        clay=5,
                        bulk_density=1.48,
                        instrument = 'TDR')

        assert arrays_are_similar(WaterEC(sample_ECW_Odarslov_top), expected_results['test_sample_ECW_Odarslov_top'])


def test_sample_ECW_Hil_ex():
        sample_ECW_Hil_ex = Soil( bulk_ec =        [0.025, 0.038, 0.065, 0.079, 0.1  ],
                        bulk_perm = [11.5,  15,    19,    22 ,   26   ],
                        clay=0,
                        bulk_density=1.8,
                        instrument = 'TDR')

        assert arrays_are_similar(WaterEC(sample_ECW_Hil_ex), expected_results['test_sample_ECW_Hil_ex'])


def test_sample_ECW1():
        sample_ECW1 = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044])

        assert arrays_are_similar(WaterEC(sample_ECW1), expected_results['test_sample_ECW1'])


def test_sample_ECW2():
        sample_ECW2 = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044],
                        temperature = np.array([15,  10,  0,   40,   15, 15]) + 273)

        assert arrays_are_similar(WaterEC(sample_ECW2), expected_results['test_sample_ECW2'])   


################################################################################################################
####################################### PREDICT WATER FROM BULK PERM ###########################################
################################################################################################################


def test_sample_WP0():
        sample_WP0 = Soil( bulk_perm =        [10,   15, 20,   25,   7,  1,  12,  22,  5,  20,  30   ], 
                                bulk_density=1.7, texture = 'Sand', solid_perm = 5, instrument = 'GPR')

        assert arrays_are_similar(Water(sample_WP0), expected_results['test_sample_WP0'])  


def test_sample_WP0b():
      sample_WP0b = Soil(water = [    0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                  bulk_perm=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ]), 
                  solid_perm = 5)

      assert arrays_are_similar(Water(sample_WP0b), expected_results['test_sample_WP0b'])  


def test_sample_WP1():
      sample_WP1 = Soil(water =[0.05, 0.11, 0.08, np.nan, np.nan, 0.07      ], 
                  bulk_perm=   [6,    11,   9,    np.nan, 8,      8.5,    12], 
                  bulk_density=1.7,
                  instrument = 'TDR')

      assert arrays_are_similar(Water(sample_WP1), expected_results['test_sample_WP1'])  


def test_sample_WP1b():
      sample_WP1b = Soil(water =     [0.05, 0.11, 0.08, np.nan, np.nan, 0.07,   np.nan, 0.2,    0.02,   np.nan ], 
                  bulk_perm=         [6,    11,   9,    np.nan, 8,      8.5,    14,     np.nan, np.nan, 1      ], 
                  bulk_density=1.7,
                  texture = 'Sand',
                  solid_perm = 5,
                  instrument = 'TDR')

      assert arrays_are_similar(Water(sample_WP1b), expected_results['test_sample_WP1b'])  


def test_sample_WP1c():
      sample_WP1c = Soil(water =    [0.20, 0.31, 0.36, 0.38, 0.05                        ], 
                        bulk_perm=  [10,   15,   20,   25,   7,   1, 12, 22, 5, 20, 30   ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        solid_perm = 5)

      assert arrays_are_similar(Water(sample_WP1c), expected_results['test_sample_WP1c'])  


def test_sample_WP3():
      sample_WP3 = Soil(water =     [0.20, 0.30, 0.35                                             ], 
                        bulk_perm = [10,   15,   20,   8.5, 8,    1,   12,   22,   5,    20,   30   ], 
                        bulk_density=1.7,
                        texture = 'Sand',
                        solid_perm = 5,
                        instrument = 'GPR')

      assert arrays_are_similar(Water(sample_WP3), expected_results['test_sample_WP3'])  


def test_sample_WP4():
      sample_WP4 = Soil(water =            [0.20, 0.30, 0.35                                             ], 
                        bulk_perm=         [10,   15,   20,  8.5,   8,   1,   12,   22,   5,    20,   30 ], 
                        bulk_density=1.7,
                        clay = 40,
                        solid_perm = 5,
                        instrument = 'TDR')

      assert arrays_are_similar(Water(sample_WP4), expected_results['test_sample_WP4'])  


def test_sample_WP5():
      sample_WP5 = Soil(water =     [0.20, 0.30, 0.35                              ], 
                        bulk_perm = [10,   15,   20,  8.5, 8, 1, 12, 22, 5, 20, 30 ], 
                        bulk_density=1.7, alpha = 0.3, frequency_perm = [150e6])

      assert arrays_are_similar(Water(sample_WP5), expected_results['test_sample_WP5'])  


def test_sample_WP7b():
      sample_WP7b = Soil(water =     [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                  bulk_perm =        [6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5   ], 
                  frequency_perm =   [50e6, 50e6, 50e6, 200e6,  200e6,  200e6,  50e6,   50e6, 50e6,   200e6 ],
                  bulk_density=1.7,
                  texture = 'Sand',
                  solid_perm = 5)

      assert arrays_are_similar(Water(sample_WP7b), expected_results['test_sample_WP7b'])  


def test_sample_WP7c():
      sample_WP7c = Soil(water =     [0.05,   0.11,   0.08,   0.11,   np.nan, np.nan, np.nan, 0.07,   np.nan, np.nan], 
                  bulk_perm =        [6,      11,     9,      np.nan, 1,      np.nan, 8,      8.5,    8.5,    8.5   ],
                  bulk_ec =          [np.nan, np.nan, np.nan, 0.002,  np.nan, np.nan, 0.003, np.nan, np.nan, np.nan],
                  frequency_perm =   [50e6,   50e6,   50e6,   200e6,  200e6,  200e6,  50e6,   50e6,   50e6,   200e6 ],
                  frequency_ec =     [50e2,   50e2,   50e2,   200e1,  200e3,  200e4,  50e3,   50e3,   50e3,   20    ],
                  bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05)
      
      assert arrays_are_similar(Water(sample_WP7c), expected_results['test_sample_WP7c'])  


def test_sample_WP8():
      sample_WP8 = Soil( bulk_perm = [10,   15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
            bulk_density=1.7,
            texture = 'Sand',
            solid_perm = 5,
            frequency_perm =         [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6])

      assert arrays_are_similar(Water(sample_WP8), expected_results['test_sample_WP8'])  


def test_sample_WP8b():
        sample_WP8b = Soil( bulk_perm =    [10,   15,    20,    25,    7,     1,    12,    20,    5,      20,    22 ], 
                        frequency_perm = [1e6,  2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,  np.nan, 100e6, 200e6],
                        bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1)

        assert arrays_are_similar(Water(sample_WP8b), expected_results['test_sample_WP8b'])  


def test_sample_WP8c():
      sample_WP8c = Soil( bulk_perm =    [10,    15,    20,    25,    7,     1,    12,    np.nan, 5,      20,    22 ], 
                        frequency_perm = [1e6,   2e6,   2.5e6, 3e6,   3.5e6, 10e6, 25e6,  25e6,   np.nan, 100e6, 200e6],
                        bulk_density=1.7, texture = 'Clay', solid_perm = 5, water_ec = 0.1)

      assert arrays_are_similar(Water(sample_WP8c), expected_results['test_sample_WP8c']) 


def test_sample_WP9():
      sample_WP9 = Soil( bulk_perm = [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
            bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_perm = 20e6)

      assert arrays_are_similar(Water(sample_WP9), expected_results['test_sample_WP9']) 


def test_sample_WP9b():
      sample_WP9b = Soil( bulk_perm = [10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ], 
            bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_perm = 1e6)

      assert arrays_are_similar(Water(sample_WP9b), expected_results['test_sample_WP9b']) 


def test_sample_WPv():
      sample_WPv = Soil( bulk_perm = [3,    8,       15,    20,    22,    7,    12,    18     ], 
                        bulk_density=1.4, texture = 'Sand', solid_perm = 5, CEC = 1.6, frequency_perm = 50e6)

      assert arrays_are_similar(Water(sample_WPv), expected_results['test_sample_WPv']) 

################################################################################################################
####################################### PREDICT WATER FROM BULK EC #############################################
################################################################################################################

def test_sample_WEC1():
      sample_WEC1 = Soil( bulk_ec = np.array([10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                bulk_density=1.7,
                texture = 'Sand',
                water_ec = 0.1)

      assert arrays_are_similar(Water(sample_WEC1), expected_results['test_sample_WEC1']) 


def test_sample_WEC1b():
      sample_WEC1b = Soil( bulk_ec=np.array(  [10,   15,   20,   25,     7,    1,     12,    22,    5,    20,    30   ])*1e-3, 
                bulk_density=1.7, texture = 'Sand',
                water_ec = np.array( [ 0.05, 0.06, 0.07, np.nan, 0.01, 0.1]))

      assert arrays_are_similar(Water(sample_WEC1b), expected_results['test_sample_WEC1b']) 


def test_sample_WEC2():
      sample_WEC2 = Soil(water = np.array([0.20, 0.31, 0.36, 0.38, 0.05                                   ]), 
                        bulk_ec= np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20,   30   ])*1e-3, 
                        bulk_density=1.7, texture = 'Sand')

      assert arrays_are_similar(Water(sample_WEC2), expected_results['test_sample_WEC2']) 


def test_sample_WEC3():
      sample_WEC3 = Soil(water = np.array([0.20, 0.31, 0.36, 0.38, 0.05                                 ]), 
                        bulk_ec=np.array([10,   15,   20,   25,   7,    1,   12,   22,   5,   20, 30   ])*1e-3, 
                        bulk_density=1.7, water_ec=0.5, texture = 'Sand')

      assert arrays_are_similar(Water(sample_WEC3), expected_results['test_sample_WEC3']) 


def test_sample_WEC4():
      sample_WEC4 = Soil(water =    [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07], 
                  bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5 ])*1e-3, 
                  water_ec = 0.05)

      assert arrays_are_similar(Water(sample_WEC4), expected_results['test_sample_WEC4']) 


def test_sample_WEC4b():
      sample_WEC4b = Soil(water = np.array([     0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                              bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', water_ec = 0.05)

      assert arrays_are_similar(Water(sample_WEC4b), expected_results['test_sample_WEC4b']) 


def test_sample_WEC5():
      sample_WEC5 = Soil( bulk_ec=np.array([          10,   15,   20,   25,   7,    1,    12,    20,  5,    20,    22 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', water_ec = 0.01, frequency_ec = 500)

      assert arrays_are_similar(Water(sample_WEC5), expected_results['test_sample_WEC5']) 


def test_sample_WEC5b():
      sample_WEC5b = Soil( bulk_ec=np.array([           10,    0,    np.nan, np.nan,   7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.1, frequency_ec = np.array([2e6]))

      assert arrays_are_similar(Water(sample_WEC5b), expected_results['test_sample_WEC5b']) 


def test_sample_WEC6():
      sample_WEC6 = Soil(water = np.array([0.20, 0.30, 0.35]), 
            bulk_ec=np.array( [10,   15,   20,   8.5,  8,    1,    12,   22,   5,    20,   30   ])*1e-3, 
            bulk_density=1.7, clay = 40, water_ec=0.4, frequency_ec=5e3)

      assert arrays_are_similar(Water(sample_WEC6), expected_results['test_sample_WEC6']) 


def test_sample_WEC6b():
      sample_WEC6b = Soil( bulk_ec=np.array([   10,   15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
                                    water =         [0.1,  0.12], bulk_density=1.7, texture = 'Sand', water_ec = 0.1,
                        frequency_ec=np.array([1,    2,     2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]))

      print('sample_WEC6b.Lw', sample_WEC6b.Lw)
      assert arrays_are_similar(Water(sample_WEC6b), expected_results['test_sample_WEC6b'])    


def test_sample_WEC6c():
      sample_WEC6c = Soil( bulk_ec=np.array([           10,    15,    20,    25,    7,     1,    12,    20,    5,    20,    22 ])*1e-3, 
            frequency_ec = np.array([1,    2,      2.5,   3,     3.5,   10,   25,    25,    50,   100,   200]),
            bulk_density=1.7, texture = 'Clay', water_ec = 0.1)

      assert arrays_are_similar(Water(sample_WEC6c), expected_results['test_sample_WEC6c'])  


def test_sample_WEC7():
      sample_WEC7 = Soil(water =           [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan], 
                        bulk_ec =np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                        bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05,
                        frequency_ec =    [50e1, 50e2, 50e2, 200e2,  200e2,  200e2,  50e2,   50e1, 50e1,   200e2] )
      
      assert arrays_are_similar(Water(sample_WEC7), expected_results['test_sample_WEC7'])  


def test_sample_WEC7b():
      sample_WEC7b = Soil(water = np.array(    [0.05, 0.11, 0.08, 0.11,   np.nan, np.nan, np.nan, 0.07, np.nan, np.nan]), 
                              bulk_ec=np.array([6,    11,   9,    np.nan, 1,      np.nan, 8,      8.5,  8.5,    8.5 ])*1e-3, 
                              bulk_density=1.7, texture = 'Sand', solid_perm = 5, water_ec = 0.05,
                        frequency_ec =         [50,   5,   50,    200e2,  200e2,  200e2,  50e2,   50e1, 50,     20])

      assert arrays_are_similar(Water(sample_WEC7b), expected_results['test_sample_WEC7b'])  


def test_sample_WECv():
      sample_WECv = Soil( bulk_ec=np.array([3,    8,    15,   20,   22,    7,     12,    18,   10,   2     ])*1e-3, 
                  frequency_ec = np.array( [50,   500,  5000, 2000, 50000, 55000, 25000, 100,  10,   50]),
                  bulk_density=1.4, texture = 'Sand', CEC = 1.6, water_ec = 0.1)

      assert arrays_are_similar(Water(sample_WECv), expected_results['test_sample_WECv']) 



