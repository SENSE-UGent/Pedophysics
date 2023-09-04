import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.insert(0, 'C:\\Users\\gmendoza\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')
sys.path.insert(0, 'C:\\Users\\mendo\\OneDrive - UGent\\Documentos\\PhD\\Pedophysics_code')

from pedophysics.simulate import Soil
from pedophysics.predict import WaterEC
from pedophysics.pedophysical_models.bulk_ec import Rhoades
from pedophysics.pedophysical_models.bulk_perm import Hilhorst


##############
########################## Testing Water EC from Bulk EC with examples from Brovelli & Cassiani 2011 ###############################
##############

print("################## Example DR_SCL ####################") 

DR_SCL = Soil( bulk_ec = [0, 1.6e-3, 4e-3, 9e-3, 1.5e-2, 2e-2],
                water =  [0, 0.076,  0.15, 0.23, 0.3,    0.38])

wec2 = WaterEC(DR_SCL) # 0.0714
print("wec2", wec2)

DR_SCL.info.to_excel('DR_SCL_info.xlsx')
DR_SCL.df.to_excel('DR_SCL_df.xlsx')

rhoadesEC =        Rhoades(DR_SCL.df.water, 0.068793, 0.0, 1.000398, 0.37895)
#----------------------------------------------------------------------------------------
fig = plt.figure()
axdr = fig.add_subplot(1, 1, 1)
aa = 0.5
ss = 50

axdr.set_title('DR_SCL')
axdr.scatter(DR_SCL.df.water, DR_SCL.df.bulk_ec, marker='D', color='black', s=ss)
axdr.plot(DR_SCL.df.water, rhoadesEC, 'ro', alpha=aa, markersize=8)
plt.show()


print("################## Example DR_L ####################") 

DR_L = Soil( bulk_ec =  [0, 7*10**-3, 1.3*10**-2, 2*10**-2, 3*10**-2, 3.3*10**-2],
                water = [0, 0.088,    0.18,       0.26,     0.35,     0.44      ])

wec2b = WaterEC(DR_L) # 0.0714
print("wec2b", wec2b)

DR_L.info.to_excel('DR_L_info.xlsx')
DR_L.df.to_excel('DR_L_df.xlsx')

rhoadesEC2b = Rhoades(DR_L.df.water, 0.093240381, 0.0027970, 0.344834, 0.62287)
#-------------------------------------------------------------------------------------
fig = plt.figure()
axdr = fig.add_subplot(1, 1, 1)

axdr.set_title('DR_L')
axdr.scatter(DR_L.df.water, DR_L.df.bulk_ec, marker='D', color='black', s=ss)
axdr.plot(DR_L.df.water, rhoadesEC2b, 'ro', alpha=aa, markersize=8)

plt.show()

print("################## Example DR_S ####################") 
                         #          0      1      2     3      4     5      6      7      8      9      10     11
DR_S = Soil( bulk_ec =  [0, 8*10**-4, 3*10**-3, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                water = [0, 0.072,      0.144,    0.22,     0.29,       0.36      ])

wec2c = WaterEC(DR_S) # 0.065
print("wec2c", wec2c)

DR_S.info.to_excel('DR_S_info.xlsx')
DR_S.df.to_excel('DR_S_df.xlsx')

rhoadesEC2c =        Rhoades(DR_S.df.water, 0.06344, 0.0, 2.19037, 0.018659)
#-----------------------------------------------------------------------------------
fig = plt.figure()
axdr = fig.add_subplot(1, 1, 1)

axdr.set_title('DR_S')
axdr.scatter(DR_S.df.water, DR_S.df.bulk_ec, marker='D', color='black', s=ss)
axdr.plot(DR_S.df.water, rhoadesEC2c, 'ro', alpha=aa, markersize=8)

plt.show()

print("################## Example DR_Sa ####################") 
                         #          0      1      2     3      4     5      6      7      8      9      10     11
DR_Sa = Soil( bulk_ec =  [0, 8*10**-4, np.nan, 6.5*10**-3, 1.3*10**-2, 1.8*10**-2],
                water =  [0, 0.072,    0.144,  np.nan,     0.29,       0.36      ])

wec2ca = WaterEC(DR_Sa) # 0.065
print("wec2ca", wec2ca)

DR_Sa.info.to_excel('DR_Sa_info.xlsx')
DR_Sa.df.to_excel('DR_Sa_df.xlsx')

rhoadesEC2ca =        Rhoades(DR_Sa.df.water, 0.6692, 0.0, 1.0, 0.38)
#-----------------------------------------------------------------------------------
fig = plt.figure()
axdr = fig.add_subplot(1, 1, 1)

axdr.set_title('DR_Sa')
axdr.scatter(DR_Sa.df.water, DR_Sa.df.bulk_ec, marker='D', color='black', s=ss)
axdr.plot(DR_Sa.df.water, rhoadesEC2ca, 'ro', alpha=aa, markersize=8)

plt.show()

print("################## Example Odarslov_top ####################") 
                         #       0     1     2     3     4     5      6      7      8      9      10     11
Odarslov_top = Soil( bulk_ec =  [0.02, 0.03, 0.04, 0.05, 0.06   ],
                    bulk_perm = [11.5, 14.8, 17,   20,   23     ],
                    clay=5,
                    bulk_density=1.48,
                    instrument = 'TDR')

wec2d = WaterEC(Odarslov_top) # 0.29
print("wec2d", wec2d)

Odarslov_top.info.to_excel('Odarslov_top_info.xlsx')
Odarslov_top.df.to_excel('Odarslov_top_df.xlsx')

HilhorstEC2d = Hilhorst(Odarslov_top.df.bulk_ec, 0.283687902, Odarslov_top.df.water_perm, 5.98)
#---------------------------------------------------------------------------------------
fig = plt.figure()
ax2d = fig.add_subplot(1, 1, 1)

ax2d.set_title('Odarslov_top')
ax2d.scatter(Odarslov_top.df.bulk_ec, Odarslov_top.df.bulk_perm, marker='D', color='black', s=ss)
ax2d.plot(Odarslov_top.df.bulk_ec, HilhorstEC2d, 'ro', alpha=aa, markersize=8)

plt.show()


print("################## Example Hil_ex ####################") 
                         #       0     1     2     3     4     5      6      7      8      9      10     11
Hil_ex = Soil( bulk_ec =        [0.025, 0.038, 0.065, 0.079, 0.1  ],
                    bulk_perm = [11.5,  15,    19,    22 ,   26   ],
                    clay=0,
                    bulk_density=1.8,
                    instrument = 'TDR')

wec2e = WaterEC(Hil_ex ) # 0.4
print("wec2e", wec2e)

Hil_ex.info.to_excel('Hil_ex_info.xlsx')
Hil_ex.df.to_excel('Hil_ex_df.xlsx')

HilhorstEC2e = Hilhorst(Hil_ex.df.bulk_ec, 0.4275173, Hil_ex.df.water_perm, 7.21040)
#-----------------------------------------------------------------------------------------
fig = plt.figure()
ax2e = fig.add_subplot(1, 1, 1)

ax2e.set_title('Hil_ex')
ax2e.scatter(Hil_ex.df.bulk_ec, Hil_ex.df.bulk_perm, marker='D', color='black', s=ss)
ax2e.plot(Hil_ex.df.bulk_ec, HilhorstEC2e, 'ro', alpha=aa, markersize=8)

plt.show()


print("################## Example Hil_exa ####################") 
                         #       0     1     2     3     4     5      6      7      8      9      10     11
Hil_exa = Soil( bulk_ec =        [0.025, 0.038, 0.065, 0.079, 0.1  ],
                    bulk_perm = [11.5,  15,    19,    22 ,   26   ],
                    clay=0,
                    bulk_density=1.8,
                    instrument = 'TDR')

wec2ea = WaterEC(Hil_exa ) # 0.4
print("wec2ea", wec2ea)

Hil_exa.info.to_excel('Hil_exa_info.xlsx')
Hil_exa.df.to_excel('Hil_exa_df.xlsx')

HilhorstEC2ea = Hilhorst(Hil_exa.df.bulk_ec, 0.4275173, Hil_exa.df.water_perm, 7.21040)
#-----------------------------------------------------------------------------------------
fig = plt.figure()
ax2e = fig.add_subplot(1, 1, 1)

ax2e.set_title('Hil_exa')
ax2e.scatter(Hil_exa.df.bulk_ec, Hil_exa.df.bulk_perm, marker='D', color='black', s=ss)
ax2e.plot(Hil_exa.df.bulk_ec, HilhorstEC2ea, 'ro', alpha=aa, markersize=8)

plt.show()


##############
########################## Testing water_ec from salinity ###############################
##############

print("################## Example1 ####################") 
                         #           0    1    2    3     4   5    6    7    8     9    10   11
sample1S = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044])

wec1S = WaterEC(sample1S)
print("wec1S", wec1S)

sample1S.info.to_excel('sample1S_info.xlsx')
sample1S.df.to_excel('sample1S_df.xlsx')

print("################## Example1b ####################") 
                         #           0    1    2    3     4   5    6    7    8     9    10   11
sample1Sb = Soil( salinity = [0.008, 0.017, 0.004, 0.026, 0.005, 0.003, 0.012, 0.021, 0.019, 0.024, 0.026, 0.044],
                temperature = np.array([15,  10,  0,   40,   15, 15]) + 273  
)
wec1Sb = WaterEC(sample1Sb)
print("wec1Sb", wec1Sb)

sample1Sb.info.to_excel('sample1Sb_info.xlsx')
sample1Sb.df.to_excel('sample1Sb_df.xlsx')
#---------------------------------------------------------------------------------------------
aa = 0.5
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.set_title('Example1 and 1b')
ax.plot(sample1S.salinity, sample1S.df.water_ec, 'bo', alpha=aa, markersize=8)
ax.plot(sample1Sb.salinity, sample1Sb.df.water_ec, 'ro', alpha=aa, markersize=8)
ax.set_ylim(0, 0.6)
ax.set_xlim(0, 0.04)

plt.show()