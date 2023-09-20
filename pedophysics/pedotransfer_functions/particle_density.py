def Schjonnen(clay, org, densorg = 1.4, denspart = 2.65, densclay = 2.86):
    """
    Calculate the soil particle density using the Schjonnen model.
    
    This function determines the particle density of soil based on the content of clay, organic matter, 
    and specific densities of organic matter, soil particles, and clay. It implements the Schjonnen model 
    (referenced year or name of the paper).
    
    Parameters
    ----------
    clay : array-like
        Soil clay content
    org : array-like
        Soil volumetric organic matter
    densorg : float, optional
        Density of organic matter [g/cm^3]. Default is 1.4 g/cm^3.
    denspart : float, optional
        Density of soil particles [g/cm^3]. Default is 2.65 g/cm^3.
    densclay : float, optional
        Density of clay particles [g/cm^3]. Default is 2.86 g/cm^3.
        
    Returns
    -------
    pd : array-like
        Soil particle density [g/cm^3].
        
    References
    ----------
    P. Schjønning, R.A. McBride, T. Keller, P.B. Obour,
    Predicting soil particle density from clay and soil organic matter contents,
    Geoderma, Volume 286, 2017, Pages 83-87, ISSN 0016-7061,
    https://doi.org/10.1016/j.geoderma.2016.10.020.
    
    Example
    -------
    >>> Schjonnen(25, 5)
    2.606
    """   
    a = 1.127
    b = 0.373
    c = 2.648
    d = 0.209

    clay = clay/100
    org = org/100

    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    return pd
