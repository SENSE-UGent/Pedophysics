import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
#sys.path.insert(0, 'C:\\Users\\gasto\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import Salinity, WaterEC


##############
########################## Testing salinity from water_ec ###############################
##############

print("################## Example1 ####################") 
                         #           0    1    2    3     4   5    6    7    8     9    10   11
sample1S = Soil( water_ec = np.array([100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3, 
)

sal1 = Salinity(sample1S) # [0.00847 0.01721 0.00419 0.02614 0.00504 0.00334 0.01281 0.02165 0.01898 0.02434 0.02614 0.04446]
print("sal1", sal1)

sample1S.info.to_excel('sample1S_info.xlsx')
sample1S.df.to_excel('sample1S_df.xlsx')

print("################## Example1b ####################") 
                         #              0    1    2    3     4   5    6    7    8     9    10   11
sample1Sb = Soil( water_ec = np.array(  [100, 200, 50,  300,  60, 40,  150, 250, 220,  280, 300, 500 ])*10**-3,
                temperature = np.array([15,  10,  0,   40,   15, 15]) + 273 
)

sal1b = Salinity(sample1Sb) # [0.0109  0.02594 0.00938 0.01969 0.00647 0.00429 0.01281 0.02165 0.01898 0.02434 0.02614 0.04446]
print("sal1b", sal1b)

sample1Sb.info.to_excel('sample1Sb_info.xlsx')
sample1Sb.df.to_excel('sample1Sb_df.xlsx')

#------------------------------------------------------

aa = 0.5
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1 and 1b')
ax.plot(sample1S.water_ec, sample1S.df.salinity, 'bo', alpha=aa, markersize=8)
ax.plot(sample1Sb.water_ec, sample1Sb.df.salinity, 'ro', alpha=aa, markersize=5)
ax.set_xlim(0, 0.6)
ax.set_ylim(0, 0.04)

plt.show()


###########################################################################

print("################## Example S ####################")     
                         #   0     1    2     3       4       5       6       7
sampleS = Soil( bulk_ec =  [0.02, 0.03, 0.04, 0.05, 0.06   ],
                bulk_perm = [11.5, 14.8, 17,   20,   23    ],
                clay=5,
                bulk_density=1.48,
                instrument = 'TDR')

print(sampleS.df)
print("predict.WaterEC(sampleS) :", WaterEC(sampleS)) 
print("predict.Salinity(sampleS) :", Salinity(sampleS)) 

sampleS.info.to_excel('sampleS_info.xlsx')
sampleS.df.to_excel('sampleS_df.xlsx')
sampleS.info.water_ec

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ss=150
ax.scatter(sampleS.df.bulk_ec, sampleS.df.perm, marker = "^", color='blue', alpha=1, s=ss)
ax.scatter(sampleS.df.bulk_ec, Hilhorst(sampleS.df.bulk_ec, sampleS.df.water_ec, sampleS.df.bulk_perm, sampleS.df.offset_perm), color='navy', alpha=1, s=ss)


fig.set_figheight(5.5)
fig.set_figwidth(8)

ax.set_ylim(0, 20)
ax.set_xlim(0, 15e-3)
ax.grid(True)
ax.set_ylabel('Water [%]', fontsize = 22)
ax.set_xlabel('Bulk real EC [S/m]', fontsize = 22)
ax.tick_params(axis='y', labelsize=16) 
ax.tick_params(axis='x', labelsize=16) 

plt.savefig('ExS1')
plt.show()