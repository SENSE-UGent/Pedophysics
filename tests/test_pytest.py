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