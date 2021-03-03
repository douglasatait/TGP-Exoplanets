#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 13 16:22:07 2021

@author: luciescharre
"""

# planet size
import batman
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.ticker as ticker
plt.rcParams["font.family"] = "serif"


params = batman.TransitParams()
params.t0 = 0.   
params.a = 8.240983371                  #semi-major axis (in units of stellar radii)
params.w =  90.                       #longitude of periastron (in degrees)
params.limb_dark = "quadratic"        #limb darkening model
params.u = [0.525,0.188]      #limb darkening coefficients [u1, u2, u3, u4]
params.inc = 84.7                      #orbital inclination (in degrees)
params.ecc = 0.121                       #eccentricity
params.rp = 0.135    
params.per = 2.152060692723822


t = np.linspace(-0.05, 0.05, 100)

m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve

plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
plt.plot(t, flux1, color='tomato')
plt.title("Literature values transit curve")
plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()

#params = batman.TransitParams()
#params.t0 = 0.                        #time of inferior conjunction
#params.per = 0.5          #orbital period days
#params.a = 8.240983371                  #semi-major axis (in units of stellar radii)
#params.w =  0.                       #longitude of periastron (in degrees)
#params.limb_dark = "nonlinear"        #limb darkening model
#params.u = []
#params.u = [0.5, 0.1, 0.1, -0.1]      #limb darkening coefficients [u1, u2, u3, u4]
#params.rp = 0.5                      #planet radius (in units of stellar radii)
   

#t = np.linspace(-0.05, 0.05, 100)


#params.inc = 90                      #orbital inclination (in degrees)
#params.ecc = 0                       #eccentricity



#params.rp = 0.1    
#params.per =  2.152060692723822           #orbital period days
#params.a = 8.240983371                  #semi-major axis (in units of stellar radii)


m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve

params.rp = 0.2    
#params.per =  2.152060692723822           #orbital period days
#params.a = 8.240983371                  #semi-major axis (in units of stellar radii)


m = batman.TransitModel(params, t)    #initializes model
flux2 = m.light_curve(params)          #calculates light curve


plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux1, label="Rp = 0.135",color='tomato')
plt.plot(t, flux2, label="Rp = 0.2")
#plt.plot(t, flux3, label="ecc: 0.5", color = "tomato")
#plt.plot(t, flux4, label="ecc: 0.8", color = "teal")
plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()


# semi major axis 

params = batman.TransitParams()
params.t0 = 0.   
params.a = 8.240983371                  #semi-major axis (in units of stellar radii)
params.w =  90.                       #longitude of periastron (in degrees)
params.limb_dark = "quadratic"        #limb darkening model
params.u = [0.525,0.188]      #limb darkening coefficients [u1, u2, u3, u4]
params.inc = 84.7                      #orbital inclination (in degrees)
params.ecc = 0.121                       #eccentricity
params.rp = 0.135    
params.per = 2.152060692723822

t = np.linspace(-0.05, 0.05, 100)


#params.inc = 90                      #orbital inclination (in degrees)

  
#params.per =  2.152060692723822           #orbital period days
#params.a = 8.240983371                  #semi-major axis (in units of stellar radii)


m = batman.TransitModel(params, t)    #initializes model
flux3 = m.light_curve(params)          #calculates light curve


#params.per =  2.15           #orbital period days
params.a = 10                  #semi-major axis (in units of stellar radii)


m = batman.TransitModel(params, t)    #initializes model
flux4 = m.light_curve(params)          #calculates light curve


params = batman.TransitParams()
params.t0 = 0.   
params.a = 7.7                 #semi-major axis (in units of stellar radii)
params.w =  90.                       #longitude of periastron (in degrees)
params.limb_dark = "quadratic"        #limb darkening model
params.u = [0.525,0.188]      #limb darkening coefficients [u1, u2, u3, u4]
params.inc = 84.7                      #orbital inclination (in degrees)
params.ecc = 0.121                       #eccentricity
params.rp = 0.135    
params.per = 2.152060692723822



plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux3, label="a = 7.7", color = "tomato")
plt.plot(t, flux4, label="a = 10")
#plt.plot(t, flux3, label="ecc: 0.5", color = "tomato")
#plt.plot(t, flux4, label="ecc: 0.8", color = "teal")
plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()


params.w = 0                      #orbital inclination (in degrees)
params.ecc = 0                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve



params.ecc = 0.121                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux2 = m.light_curve(params)          #calculates light curve


params.ecc = 0.5                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux3 = m.light_curve(params)          #calculates light curve

params.ecc = 0.3                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux4 = m.light_curve(params)          #calculates light curve


plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux1, label="ecc = 0, circular orbit", color = "green")
plt.plot(t, flux2, label="ecc = 0.121",color='tomato')
plt.plot(t, flux4, label="ecc: 0.3")
plt.plot(t, flux3, label="ecc = 0.5")

#plt.title("Longitude of periastron "+ '\u03C9'+'\u0305'+ "= 0")
plt.vlines(0,1,0.982)
plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()


params.w = 90                      #orbital inclination (in degrees)
params.ecc = 0                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve



params.ecc = 0.121                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux2 = m.light_curve(params)          #calculates light curve


params.ecc = 0.5                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux3 = m.light_curve(params)          #calculates light curve

params.ecc = 0.3                       #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux4 = m.light_curve(params)          #calculates light curve

plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux1, label="e = 0, circular orbit", color = "green")
plt.plot(t, flux2, label="e = 0.121",color='tomato')
plt.plot(t, flux4, label="e = 0.3")
plt.plot(t, flux3, label="e = 0.5")

#plt.plot(t, flux4, label="ecc: 0.8", color = "teal")
#plt.title("Longitude of periastron "+ '\u03C9'+'\u0305'+ "= 90")
plt.vlines(0,1,0.982)
plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()


params = batman.TransitParams()
params.t0 = 0.   
params.a = 8.240983371                  #semi-major axis (in units of stellar radii)
params.w =  90.                       #longitude of periastron (in degrees)
params.limb_dark = "quadratic"        #limb darkening model
params.u = [0.525,0.188]      #limb darkening coefficients [u1, u2, u3, u4]
params.inc = 84.7                      #orbital inclination (in degrees)
params.ecc = 0.121                       #eccentricity
params.rp = 0.135    
params.per = 2.152060692723822

#params.w = 90                      #orbital inclination (in degrees)


m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve



params.inc = 90                      #eccentricity


m = batman.TransitModel(params, t)    #initializes model
flux2 = m.light_curve(params)          #calculates light curve



plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux1, label="i = 84.7°", color = "tomato")
plt.plot(t, flux2, label="i = 90°")
#plt.plot(t, flux4, label="ecc: 0.8", color = "teal")

plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()




params = batman.TransitParams()
params.t0 = 0.   
params.a = 8.240983371                  #semi-major axis (in units of stellar radii)
params.w =  90.                       #longitude of periastron (in degrees)
params.limb_dark = "quadratic"        #limb darkening model
params.u = [0.525,0.188]      #limb darkening coefficients [u1, u2, u3, u4]
params.inc = 84.7                      #orbital inclination (in degrees)
params.ecc = 0.121                       #eccentricity
params.rp = 0.135    
params.per = 2.152060692723822

#params.w = 90                      #orbital inclination (in degrees)


m = batman.TransitModel(params, t)    #initializes model
flux1 = m.light_curve(params)          #calculates light curve



params.limb_dark = "linear"        #limb darkening model
params.u = [0.525]      #limb darkening coefficients [u1, u2, u3, u4]


m = batman.TransitModel(params, t)    #initializes model
flux2 = m.light_curve(params)          #calculates light curve

params.limb_dark = "uniform"        #limb darkening model
params.u = []      #limb darkening coefficients [u1, u2, u3, u4]


m = batman.TransitModel(params, t)    #initializes model
flux3 = m.light_curve(params)          #calculates light curve

plt.minorticks_on()
plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)

plt.plot(t, flux3, label="uniform")
plt.plot(t, flux2, label="linear")
plt.plot(t, flux1, label="quadratic")


plt.xlabel("Time since transit [days]")
plt.ylabel("Relative flux")
#plt.title()
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()






