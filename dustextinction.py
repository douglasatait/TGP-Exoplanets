import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
import os
import pandas
from dust_extinction.averages import (GCC09_MWAvg, RL85_MWGC, I05_MWAvg, CT06_MWLoc, CT06_MWGC,F11_MWGC)

true_col_data = open("intrinsic_colours.txt", "r")
B_V =[]
U_B=[]
V_I =[]
V_R = []
sp_typ = []

for line in true_col_data.readlines() : 
   if line.startswith("#") : 
      pass
   else : 
      token = line.split(',')
      sp_typ.append(token[0]) 
      B_V.append(float(token[2])) 
      U_B.append(float(token[1])) 
      V_R.append(float(token[3]))
      V_I.append(float(token[4]))  


models =[GCC09_MWAvg]
fig, ax = plt.subplots()
x = np.arange(0.51,0.7275,0.001)/u.micron 
#the start and stop values are the wavelengths covered by the colour passbands you want to predict the spectral type using  
#see filter sheet in teams for the lower and upper limits of each band

for cmodel in models:
   ext_model = cmodel()
   indxs, = np.where(np.logical_and(
      x.value >= ext_model.x_range[0],
      x.value <= ext_model.x_range[1]))
yval = ext_model(x[indxs])
ax.plot(x[indxs], yval, label=ext_model.__class__.__name__)

yvals =yval
ax.set_yscale("log")
ax.set_xlabel(r'$\lambda$ [$\mu m$]')
ax.set_ylabel(r'$A(\lambda)/A(V)$')
ax.set_title('Dust Extinction')

#ax.legend(loc='best')
#plt.tight_layout()
#plt.show()

"""u_g1= 1.14+0.03
g_r1 = 0.29+0.03
r_i1 =0.02
u_g2 =1.35-0.04
g_r2=0.17-0.046
r_i2 =0.01
U_B_con = 0.77*u_g2-0.88
B_V_con = 0.98*g_r2+0.22
V_R_con = 1.09*r_i2+0.22
#print(U_B_con)
#print(B_V_con)
print(V_R_con)"""


#function to predict spectral type from U-B colour
"""for i in range(len(U_B)):
   for j in range (len(yvals)):
      if U_B[i]+yvals[j] >= (0.476-0.1) and U_B[i]+yvals[j] <= (0.476+0.1): #input the observed U-B colour and its error in this line 
         print("The U_B spectral type is:" ,sp_typ[i])
         print(U_B[i]+ yvals[j])
         print(U_B[i])
         print(yvals[j])"""
#function to predict spectral type from B-V colour
"""for i in range(len(B_V)):
   for j in range (len(yvals)):
      if B_V[i]+yvals[j] >= (0.636-0.3) and B_V[i]+yvals[j] <= (0.636+0.3): #input the observed B-V colour and its error in this line 
         print("The B-V spectral type is:" ,sp_typ[i])
         print(B_V[i]+ yvals[j])
         print(B_V[i])
         print(yvals[j])"""
#function to predict spectral type from V-I colour
"""for i in range(len(V_I)):
   for j in range (len(yvals)):
      if V_I[i]+ yvals[j] >= (0.759-0.2) and V_I[i]+yvals[j] <= (0.759+0.2): #input the observed V-I colour and its error in this line 
         print("The V-I spectral type is:" , sp_typ[i])
         print(V_I[i]+ yvals[j])
         print(V_I[i])
         print(yvals[j])"""

#function to predict spectral type from V-R colour
for i in range(len(V_R)):
   for j in range (len(yvals)):
      if V_R[i]+yvals[j] >= (0.390-0.25) and V_R[i]+yvals[j] <= (0.390+0.25): #input the observed V-R colour and its error in this line 
         print("The V-R spectral type is:" ,sp_typ[i])
         print(V_R[i]+ yvals[j])
         print(V_R[i])
         print(yvals[j])


    