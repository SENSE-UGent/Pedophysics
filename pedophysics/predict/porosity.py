import numpy as np
from pedophysics.predict.particle_density import ParticleDensity

def Porosity(soil):
    """
    Calculate missing values of soil.df.porosity and return

    Calculates the porosity of soil based on its bulk density and particle density. 
    If the porosity value is missing for any soil state, it is computed using the formula `1 - (bulk_density / particle_density)`.
    This function also updates the soil object's porosity information, marking it as "Calculated based on bulk density" wherever applicable.

    Parameters
    ----------
    soil : Soil Object
        An object representing the soil, which must have the following attributes:
        - df: DataFrame
            A pandas DataFrame containing the soil states with columns for `bulk_density`, `particle_density`, and `porosity`.
        - n_states: int
            The number of soil states represented in the `df`.
        - info: dict
            A dictionary containing information about the soil properties, including `porosity`.
        - roundn: int
            The number of decimal places to round the calculated porosity values to.

    Returns
    -------
    numpy.ndarray
        soil.df.porosity.values: an array of the porosity values for each soil state.

    External functions
    --------
    ParticleDensity : Calculate or set missing values of soil.df.particle_density and return

    Example
    -------
    # Assuming `soil` is a pre-defined Soil object with the required attributes
    porosity_values = Porosity(soil) 
    print(porosity_values)
    """
    ParticleDensity(soil)

    # Check for missing values
    missing_porosity_before = soil.df['porosity'].isna()

    # Calculate missing porosity values where possible
    soil.df['porosity'] = [
        round(1 - soil.df.bulk_density[x] / soil.df.particle_density[x], soil.roundn)
        if np.isnan(soil.df.porosity[x]) and not np.isnan(soil.df.bulk_density[x]) and not np.isnan(soil.df.particle_density[x])
        else soil.df.porosity[x] 
        for x in range(soil.n_states)
    ]
    missing_porosity_after = soil.df['porosity'].isna()

    # Update info for calculated porosity
    soil.info['porosity'] = [str(soil.info.porosity[x]) + (
            "--> Calculated based on bulk density"
            if missing_porosity_before[x] and not missing_porosity_after[x]
            else "--> Provide porosity or bulk_density"
            if missing_porosity_before[x] and missing_porosity_after[x]
            else "")
        if missing_porosity_before[x]
        else soil.info.porosity[x]
        for x in range(soil.n_states)]
    
    return soil.df.porosity.values


