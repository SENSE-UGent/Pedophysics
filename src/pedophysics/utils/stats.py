import numpy as np

def R2_score(actual, predicted):
    """
    Calculate the coefficient of determination (R^2) of a prediction.

    The R^2 score function computes the coefficient of determination, often used 
    to evaluate the performance of a regression model. The best possible score is 1.0.
    This function is designed to handle arrays with NaN values by ignoring such entries.

    Parameters
    ----------
    actual : array-like of shape (n_samples,)
        Ground truth (correct) target values.
    predicted : array-like of shape (n_samples,)
        Estimated targets as returned by a classifier.

    Returns
    -------
    float
        R^2 of the prediction.

    Notes
    -----
    This function works with arrays that include NaN values, ignoring such entries 
    during the computation. Therefore, 'actual' and 'predicted' arrays can have 
    missing values, but they must be of the same shape.

    Example
    -------
    >>> actual = np.array([3, -0.5, 2, 7, 4.2])
    >>> predicted = np.array([2.5, 0.0, 2.1, 7.8, 5.3])
    >>> R2_score(actual, predicted)
    0.9228556485355649
    """
    valids = ~np.isnan(actual) & ~np.isnan(predicted)
    # Calculate the total sum of squares
    ss_tot = np.sum((actual[valids] - np.mean(actual[valids])) ** 2)
    # Calculate the residual sum of squares
    ss_res = np.sum((actual[valids] - predicted[valids]) ** 2)

    # Calculate the R2 score
    r2 = 1 - (ss_res / ss_tot)
    return r2