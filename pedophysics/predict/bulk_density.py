import numpy as np
import pedotransfer_functions as ptf



def BulkDensity(soil):
### Bulk density prediction ####

    if (soil.orgm != np.nan) & (soil.OC == np.nan):
        print("OC predicted by 0.58*orgm")

        soil.set_OC(ptf.orgm2oc(soil.orgm))


    if (soil.region == "Europe") & (soil.sand == np.nan):
        print("For soils in Europe, sand and clay values are needed")

    elif (soil.region == "Europe") & (soil.clay == np.nan):
        print("For soils in Europe, sand and clay values are needed")

    elif (soil.OC != np.nan) & (soil.sand != np.nan) & (soil.clay != np.nan) & (soil.land == "cultivated topsoil") & (soil.region == "Europe"):
        print("Bulk density predicted using pedotransfer function of Hollis et al., 2012. R2 = 0.62, RMSE = 0.13. n = 333 (Europe)")

        return ptf.hollis_ct(soil.OC, soil.sand, soil.clay)

    elif (soil.OC != np.nan) & (soil.sand != np.nan) & (soil.clay != np.nan) & (soil.region == "Europe"):
        print("Bulk density predicted using pedotransfer function of Hollis et al., 2012. R2 = 0.62, RMSE = 0.15. n = 925 (Europe)")

        return ptf.hollis_mh(soil.OC, soil.sand, soil.clay)

    elif (soil.LOI != np.nan) & (soil.land == 'forest') & (soil.region == "Europe"):
        print("Bulk density predicted using pedotransfer function of De Vos et al., 2005. R2 = 0.57")

        return ptf.de_vos(soil.LOI)

    elif soil.OC != np.nan:
        print("Bulk density predicted using pedotransfer function of Abdelbaki 2018. RMSE = 0.14, R2 = 0.68, n = 45195 (USA)")

        return ptf.abdelbaki(soil.OC)

    else:
        print("To predict bulk density, please provide OC or orgm. Try: 'your_soil.OC = ...'. For more accuracy, specify land use and region.")