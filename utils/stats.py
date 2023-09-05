import numpy as np

def R2_score(actual, predicted):
    """
    """
    # Calculate the total sum of squares
    ss_tot = np.sum((actual - np.mean(actual)) ** 2)

    # Calculate the residual sum of squares
    ss_res = np.sum((actual - predicted) ** 2)

    # Calculate the R2 score
    r2 = 1 - (ss_res / ss_tot)
    
    return r2