"""
Original file created by Adam Dempsey - https://adamdempsey90.github.io/python/radial_velocities/radial_velocities.html

Edited by Douglas Tait for our TGP radial velocity fitting. 

This file is the brains of the program. It creates a star object, with its mass and RV data.
The file then contains all of the functions neccessary to calculate and return planetary parameters and create relevant graphs.

Last Edit - 03/03/2021
"""


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import scipy.signal as sig
from scipy.optimize import curve_fit,fsolve
from copy import copy
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import matplotlib.ticker as ticker

G = 4*np.pi**2 /(365.25)**2 # G in AU^3/(day^2 Solar Mass)
Gconst = 6.67e-11
Msun = 1.989e30
Mjup = 1.898e27
transit_time = 2453991.5146 - 2450000
au = 1.496597870e11


class Star():
    def __init__(self,name,mass,rv_data,mass_err,jitter):
        self.name = name
        self.mass = mass
        self.t = rv_data[:,0] - 2450000
        self.vr = rv_data[:,1]
        self.vr_err = rv_data[:,2]
        self.mass_err = mass_err
        self.vr_err = np.sqrt(self.vr_err**2 + jitter**2)
        self.jitter = jitter
        
    def plot(self):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.errorbar(self.t,self.vr,yerr=self.vr_err,fmt='o')
        ax.set_xlabel('Time  (Julian Days)',fontsize=15)
        ax.set_ylabel('Radial Velocity (m/s)',fontsize=20)
        ax.set_title(self.name + ',    M = %.2f Msun' % self.mass)
        ax.minorticks_on()
        ax.grid()
        plt.show()
    

    def get_best_n(self):
        periods = 10**np.linspace(-3,4,1000)
        omega = 2*np.pi/periods
        res=sig.lombscargle(self.t-self.t[0],self.vr,omega)
        max_pers = np.argsort(res)
        return omega[max_pers[-1]]
    
    
    def periodogram(self):
        periods = 10**np.linspace(-3,4,1000)
        omega = 2*np.pi/periods
        power=sig.lombscargle(self.t-self.t[0],self.vr,omega)
        self.periodogram_plot(periods,power)


    def periodogram_plot(self,periods,power,pmax=None,ax=None):
        if ax == None:
            fig=plt.figure()
            ax=fig.add_subplot(111)
        ax.semilogx(periods,power)
        ax.set_xlabel('Period (days)',fontsize=15)
        ax.set_ylabel('Power',fontsize=20)
        if pmax != None:
            ax.axvline(pmax,color='k',linestyle='--',linewidth=3)
            
            
    def recover_params(self,*popt):
        """Recover the mass of the planet"""

        n = popt[0]
        k = popt[2]
        w = popt[3]
        e = popt[4]

        
        if e<0:
            w += np.pi
            e *= -1
        if k<0:
            k*=-1
            w += np.pi
            
        p = 2*np.pi/n
        k1 = ((2 * np.pi * Gconst)/(p * 24 * 3600))**(1./3) * (1 - e**2)**(-1./2)
        mass_f = k/k1

        # Now make assumption that most of the mass is in the star
        mp = mass_f * (self.mass * Msun)**(2./3) 
        mp /= Mjup
        a = ((((p*24*3600)**2)*Gconst*self.mass*Msun)/(4*np.pi**2))**(1/3)
        a = a/au
        # a = (G*self.mass/n**2)**(1./3)
        return mp,e,p,w,a
        
    
    def get_uncertainties(self,pcov,*popt):
        mp,e,p,w,a = self.recover_params(*popt)
        
        ms = self.mass
        ms_err = self.mass_err
       
        n = popt[0]
        k = popt[2]
        
        n_err = np.sqrt(pcov[0,0])
        k_err = np.sqrt(pcov[2,2])
        w_err = np.sqrt(pcov[3,3])
        e_err = np.sqrt(pcov[4,4])

        p_err = p*np.abs(n_err/n)
        a_err = a*np.sqrt((2*n_err/(3*n))**2 + (ms_err/(3*ms))**2)
        mp_err = mp*np.sqrt( (2*ms_err/(3*ms))**2 + (k_err/k)**2 + (n_err/(3*n))**2 + (e*e/(1-e*e))**2 *(e_err/e)**2)
        
        return mp_err,e_err,p_err,w_err,a_err
    
    
    def output_results(self,mp,mp_err,e,e_err,p,p_err,a,a_err):
    
        outstr = 'Mp > {:.3f} +\- {:.4f} MJ\n'.format(mp,mp_err)
        outstr += 'e = {:.3f} +\- {:.4f}\n'.format(e,e_err)
        outstr += 'P =  {:.3f} +\- {:.4f} days\n'.format(p,p_err)
        outstr += 'a = {:.3f} +\- {:.4f} AU\n'.format(a,a_err)
    
        print(self.name)
        print(outstr)



def solve_kep_eqn(l,e):
    """ Solve Keplers equation x - e*sin(x) = l for x"""
    try:
        l[0]
        res = np.zeros(l.shape)
        for i,li in enumerate(l):
            tmp,= fsolve(lambda x: x-e*np.sin(x) - li,li)
            res[i] = tmp
    except IndexError:
        res, = fsolve(lambda x: x - e*np.sin(x)-l,l)

    return res
        
    
        
def load_single_star(fname):
    with open(fname,'r') as f:
        header = f.readline()
    
    rv_data = np.loadtxt(fname)
    line = header.strip().split('#')[-1].strip()
    line = line.split('Mass')
    name = line[0]
    line = line[-1].split('Jitter')
    mass = line[0]
    jitter = float(line[-1])
    mass = float(line[0].split('(')[0])
    mass_err = float(line[0].split('(')[-1].split(')')[0])
    return Star(name,mass,rv_data,mass_err=mass_err,jitter=jitter)
    
    
def plotting(star, xdata1, ydata1, xdata2, ydata2, photophase):
    
    t = star.t
    vr = star.vr 
    vr_err = star.vr_err
    
    plt.rcParams["font.family"] = "serif"
    fig = plt.figure(figsize=(20,14))
    gs = gridspec.GridSpec(3,4)
    
    ax1 = fig.add_subplot(gs[:2,:])
    ax2 = fig.add_subplot(gs[-1,:])
    
    ax1.set_title('WASP-2: Radial Velocity Curve Fit', fontsize=35, pad=20)
    ax1.xaxis.set_minor_locator(AutoMinorLocator())
    ax1.yaxis.set_minor_locator(AutoMinorLocator())
    ax1.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True, labelsize=20)
    ax1.tick_params(which='minor', direction='in', length=4, top=True, right=True)

    ax2.xaxis.set_minor_locator(AutoMinorLocator())
    ax2.yaxis.set_minor_locator(AutoMinorLocator())
    ax2.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True, labelsize=20)
    ax2.tick_params(which='minor', direction='in', length=4, top=True, right=True)

    ax1.errorbar(t, vr, yerr=vr_err, fmt='o', ms=7.5, ecolor='black', marker='o', mew=1.5, mec='blue', mfc='blue')
    ax1.plot(xdata1, ydata1, label='Fitted Curve', color='darkorange')
    ax1.set_xlabel(r"Barycentric JD$ - 2450000.0$", fontsize=25)
    ax1.set_ylabel(r"Radial Velocity " +'\n' r"$[ms^{-1}]$", fontsize=25)

    ax2.errorbar(photophase, vr, yerr=vr_err, fmt='o', ms=7.5, ecolor='black', marker='o', mew=1.5, mec='blue', mfc='blue')
    ax2.plot(xdata2, ydata2, label='Fitted Curve', color='darkorange')
    ax2.set_xlabel("Photometric Phase", fontsize=25)
    ax2.set_ylabel(r"Radial Velocity " + '\n' r"$[ms^{-1}]$", fontsize=25)
    fig.tight_layout(pad=3.0)
    # plt.savefig('rad_vel_fit_phase.png')
    plt.show()
    
    
def plot_residuals(star,residuals):

    t = star.t
    vr = star.vr 
    vr_err = star.vr_err
             
    plt.rcParams["font.family"] = "serif"
    fig = plt.figure(figsize=(20,10))
    gs = gridspec.GridSpec(3,4)
    ax = fig.add_subplot(gs[:2,:])
    plt.subplots_adjust(hspace=0)
    fig.suptitle('WASP-2: Radial Velocity Curve Fit Residuals', fontsize=40)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True, labelsize=20)
    ax.tick_params(which='minor', direction='in', length=4, top=True, right=True)

    ax.errorbar(t,residuals,yerr=vr_err, fmt='o', ms=7.5, label='Radial Velocity Data', ecolor='black', marker='o', mew=1.5, mec='blue', mfc='blue')

    ax.set_xlabel(r"Barycentric JD$ - 2450000.0$", fontsize=25)
    ax.set_ylabel(f"Radial Velocity \n Residuals", fontsize=25)

    ax.axhline(0,color='k',linewidth=1)
    # plt.savefig('rad_vel_residuals.png')
    plt.show()
    