"""
    ...

    :AUTHOR: Gaston Mendoza Veirana
    :CONTACT: gaston.mendozaveirana@ugent.be

    :REQUIRES: numpy, scipy
"""
def Schjonnen(clay, org, densorg = 1.4, denspart = 2.65, densclay = 2.86):
    """
        Schjonnen et al. 2017

        Parameters
        ----------

        clay: float
            Soil volumetric clay content [m**3/m**3]

        orgm: float
            Soil volumetric organic matter [m**3/m**3]

        Returns
        -------
        Soil particle density [kg/m**3]
    """       
    a = 1.127
    b = 0.373
    c = 2.648
    d = 0.209

    somr = (org*densorg)/(org*densorg + (1-org)*denspart)
    claymass = (clay*densclay)/(clay*densclay + (1-clay)*denspart)
    pd = ((somr/(a+b*somr)) + (1-somr)/(c+d*claymass))**-1
    return pd
