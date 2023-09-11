import numpy as np
import sys
sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.predict import BulkEC, BulkPerm
from pedophysics.simulate import Soil

from pedophysics.pedophysical_models.bulk_ec import Rhoades


def arrays_are_similar(a, b):
    # Check if the two arrays are of the same shape
    if a.shape != b.shape:
        return False
    
    # Check if the non-NaN elements are close to each other
    non_nan_match = np.isclose(a[~np.isnan(a)], b[~np.isnan(b)]).all()

    # Check if the NaN locations are the same in both arrays
    nan_match = np.isnan(a) == np.isnan(b)

    return non_nan_match and nan_match.all()

##############
########################## Testing DC frequency, fitting and non-fitting ###############################
##############

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


        