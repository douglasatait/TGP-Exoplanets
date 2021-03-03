"""
Original file created by Adam Dempsey - https://adamdempsey90.github.io/python/radial_velocities/radial_velocities.html

Edited by Douglas Tait for our TGP radial velocity fitting. 

This file reads in the data and performs the radial velocity fit using scipy.optimize.cure_fit().
It then prints out a list of the fitted and calculated planetary properties.
It also produces graphs of the curve fit, photophase fit and residuals.

Last Edit - 03/03/2021
"""


import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit,fsolve
from rv import *


au = 1.496e11
rwasp = 0.78
# rwasp = 0.877
rsun = 6.960e8
transit_time = 2453991.5146 - 2450000
systemer = 10
# transit_time = 2453991.5146

def fitting_function(t,n,tau,k,w,e):
    """Obtain the radial velocity due to a single planet.
        t = time of measurement,
        n = angular frequency of planet,
        tau = time of pericenter passage,
        k = amplitude of radial velocity (depends on planet mass and eccentricity),
        w = related to the argument of pericenter by a shift of pi.
        e = eccentricity of orbit
        The radial velocity at time t is given by
        vr = k*(cos(f + w)+e*cos(w)),
        where f is related to the number of periods since pericenter passage, n*(t-tau)
    """



    e_anom = solve_kep_eqn(n*(t-tau),e)
    f = 2*np.arctan2(np.sqrt(1+e)*np.sin(e_anom*.5),np.sqrt(1-e)*np.cos(e_anom*.5))

    return k*(np.cos(f + w) + e*np.cos(w))



def fit_data(star,fitting_function):
    """ Fit the RV data with the given fitting function """

    # Set initial guesses to the parameters
    n0 = star.get_best_n()    
    k0 = np.max(star.vr)
    tau0 = star.t[star.vr == k0]
    if len(tau0) == 0:
        tau0 = star.t[0]
    w0 = np.pi * 0
    e0 = 0    
    initial_guess1 = [n0,tau0[0],k0,w0,e0]


    # Now do the fit
    # Pass the time and radial velocities to curve_fit 
    popt,pcov = curve_fit(fitting_function, star.t, star.vr, sigma=star.vr_err,absolute_sigma=True, p0=initial_guess1)
    perr = np.sqrt(np.diag(pcov))
    mp,e,p,w,a= star.recover_params(*popt)
    mp_err,e_err,p_err,w_err,a_err = star.get_uncertainties(pcov,*popt)

    t_fit = np.linspace(int(star.t[0]),int(star.t[-1]) + 1,int(1e3))
    vr_fit = np.array([fitting_function(x,*popt) for x in t_fit])
    residuals = star.vr - np.array([fitting_function(x,*popt) for x in star.t])
    
    fparam = np.array(['n','tau','k','w','e'])
    cparamname = np.array(['Mp','e','p','w','a'])
    cparam = np.array([mp,e,p,w,a])
    cparamerr = np.array([mp_err,e_err,p_err,w_err,a_err])
        
    photophase = ((star.t - transit_time/4) % (p)) - (p/2) 
    t_fit2 = np.linspace(int(photophase[0])-1,int(photophase[-1]) + 1,int(1e3))
    vr_fit2 = np.array([fitting_function(x,*popt) for x in t_fit2])

    chisq = np.sum(residuals**2/star.vr_err**2)
    rchisq = chisq/(len(star.vr)-5)

    # Output the results by calling the star.output_results
    for i in range(len(popt)):
        print(f'Fitted paramter, {fparam[i]} = {popt[i]} +/- {perr[i]}\n')
    print('\n')
    for i in range(len(cparam)):
        print(f'Calculated paramter, {cparamname[i]} = {cparam[i]} +/- {cparamerr[i]}\n')
    print('\n')
    print(f'Chi Squared: {chisq} \nReduced Chi Squared: {rchisq}\n')
    
    # Finally plot the results by calling the plot data function
    plotting(star, t_fit, vr_fit, t_fit2, vr_fit2, photophase)
    plot_residuals(star,residuals)
    
    return mp,e,p,w,a,mp_err,e_err,p_err,w_err,a_err


plt.close('all')
# Load the data file in here

star = load_single_star('radveldata.txt')
# star = load_single_star('radveldata2.txt')
# star = load_single_star('hd10442.dat')


# Now fit the data with the fit_data function defined above
mp,e,p,w,a,mp_err,e_err,p_err,w_err,a_err = fit_data(star, fitting_function)

awasp = (a * au) / (rwasp * rsun)
awasp_err = (a_err * au) / (rwasp * rsun)
print(f'Semi-major axis / WASP2 Radius = {awasp} +/- {awasp_err}')


#FIX ME!

plt.show()


# Now repeat these steps for the other stars.
