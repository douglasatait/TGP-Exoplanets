#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 6 15:39:09 2021
Last Edit: Tues Mar 2 13:29:54 2021

@author: luciescharre/euannewlands

Description: This code reads in time-scaled catalogs from WASP-2 (WASAP-2b transit) observation. 
WASP-2 and 3 comparison field stars are used read in. Atmospheric effects are corrected for by using the field stars
where the final corrected WASP-2 flux is plotted and a light curve is fitted using the batman package.
The best fit light curve is determined and the model returns best fit parameters, used to determine the properties of WASP-2b
"""

import os
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import numpy as np
import batman
import scipy.optimize as spicy
from scipy.stats import chisquare

def chsq(oflux,eflux,ferr,params):	# arrays of both
	chisq=np.sum(((oflux-eflux)**2)/ferr**2)
	#pchisq=np.sum(((oflux-eflux)**2)/eflux**2)
	#print(chisq,len(oflux))
	#print(pchisq)
	rchisq=chisq/(len(oflux)-len(params))
	print("reduced chi-square = %s"%rchisq)
	return rchisq
	
def rJup(rstar,rerr):
	conv=0.10045		# jupiter radii in solar radii
	R_coll=0.79			# wasp2 radii in solar radii (Collier et al. 2006) 
	R_add=0.866			# wasp2 radii in solar radii (Addison et al. 2019)
	rj=rstar*R_coll/conv
	
	rjerr=rj*(rerr/rstar)
	
	return rj,rjerr

def div_err(v1,v2,e1,e2):
	# v1 is the numerator, v2 is the denominator, e1/e2 are corresponding errors
	pe1=e1/v1
	pe2=e2/v2
	p_err=np.sqrt(pe1**2 + pe2**2)	# adding percentage errors in quadrature
	err=(v1/v2)*p_err
	#print(e1[5],v1[5],pe1[5],pe2[5],p_err[5],v1[5]/v2[5],err[5])	
	return err
	
def w_mean(c1,c2,c3,e1,e2,e3):			# weighted mean of flux fluctuations from independent measurements
	vals=[c1,c2,c3]
	es=[e1,e2,e3]
	ws=[0,0,0]
	s1=[0,0,0]
	for i in range(3):
		ws[i]=1/(es[i]**2)				# smaller the error, the greater the weight
		s1[i]=vals[i]*ws[i]				# multiplying value by weight
	S1=sum(s1)
	S2=sum(ws)
	
	mean=S1/S2							# values with greater weight count more 
	merr=np.sqrt(1/S2)
	
	#print(c1[0],e1[0],mean[0],merr[0])
	
	return mean,merr
	
def scaleError(ferr,scale):
	s_err=ferr*(1/scale)			# an array of errors scaled by the same factor
	return s_err

def etest(of,ef,err):
	count=0
	for i in range(len(of)):
		if ((of[i]+err[i]>=ef[i]) and of[i]<=ef[i]) or (of[i]-err[i]<=ef[i] and of[i]>=ef[i]):
			count+=1
	frac=count/len(of)
	if frac<(2/3):			# if fraction if less than 2/3 increase error bar size
		up=True
	else:
		up=False
	
	return frac,up

def scale_err(nWFlux,MFlux,errs):
	# find the optimum scale to find 2/3 of error bars to cross
	scalar=4.5
	s_err=scaleError(errs,scalar)
	no=0.
	while no < 0.661 or no > 0.661:	# finds scale for error bars such that 2/3 cross model curve
		print('Scalar value is %s' %scalar)
		s_err=scaleError(errs,scalar)
		no,up=etest(nWFlux,MFlux,s_err)
		print(no)
		#print(no,up)
		if no > 0.661 and no < 0.671:
			print('Final scalar value is %s' %scalar)
			break
		if up==True:
			scalar-=0.001
		if up==False:
			scalar+=0.001
	return s_err
	

def read():

	NUMBER =[]
	FLUX_APER = []
	FLUXERR_APER = []
	MAG_APER = [] 
	MAGERR_APER = [] 
	X_IMAGE = []
	Y_IMAGE = []

	NUMBER_C1 =[]
	FLUX_APER_C1  = []
	FLUXERR_APER_C1  = []
	MAG_APER_C1  = [] 
	MAGERR_APER_C1  = [] 
	X_IMAGE_C1  = []
	Y_IMAGE_C1  = []

	NUMBER_C2 =[]
	FLUX_APER_C2 = []
	FLUXERR_APER_C2 = []
	MAG_APER_C2 = [] 
	MAGERR_APER_C2 = [] 
	X_IMAGE_C2 = []
	Y_IMAGE_C2 = []

	NUMBER_C3  =[]
	FLUX_APER_C3 = []
	FLUXERR_APER_C3 = []
	MAG_APER_C3 = [] 
	MAGERR_APER_C3 = [] 
	X_IMAGE_C3 = []
	Y_IMAGE_C3 = []

	NUMBER_C4  =[]
	FLUX_APER_C4 = []
	FLUXERR_APER_C4 = []
	MAG_APER_C4 = [] 
	MAGERR_APER_C4 = [] 
	X_IMAGE_C4 = []
	Y_IMAGE_C4 = []

	# read in catalog
	cata='w2b'
	#cata=str(input('catalog: '))
	# catalog choices depending on data reduction
	if cata == '1':
		catadata='catalogs'
		tfile='WASP.txt'
	elif cata=='2':
		catadata='catalogs2'
		tfile='WASP.txt'
	elif cata=='3':
		catadata='catalogs3'
		tfile='WASP.txt'
	elif cata=='lit':
		catadata='catalogsL'
		tfile='WASP_Maybe.txt'
	elif cata=='GV':
		catadata='catalogsGV'
		tfile='WASP_Maybe.txt'
	elif cata=='mb':
		catadata='catalogs3_MAYBE'
		tfile='WASP_Maybe.txt'
	elif cata=='GVf':
		catadata='catalogs_20pix_GV'
		tfile='WASP_Maybe.txt'
	elif cata=='Lf':
		catadata='catalogs_20pix_L'
		tfile='WASP_Maybe.txt'
	elif cata=='w2b':
		catadata='wasp2b_transit'
		tfile='time.txt'
	
	
	entries = os.listdir('/Users/Euan/Desktop/LC_LS/%s'%catadata)  #adjust this to your path
	entries.sort()
	for filename in entries:
		if filename.endswith('.cat'):
			with open(os.path.join('/Users/Euan/Desktop/LC_LS/%s'%catadata, filename)) as filein:
				#content = filein.read()
	
				for line in filein.readlines():                      #iterates through the lines in the file
				
					if not line.startswith('#'): 
					
						tokens = line.split()                   # breaks line into tokens seperated by ,
					
						#WASP 2
						#after telescope flip
						 #before telescope flip
						if 856<=float(tokens[5])<=861:
							if 1036<=float(tokens[6])<=1040:
								NUMBER.append(float(tokens[0]))
								FLUX_APER.append(float(tokens[1]))
								FLUXERR_APER.append(float(tokens[2]))
								MAG_APER.append(float(tokens[3]))
								MAGERR_APER.append(float(tokens[4]))
								X_IMAGE.append(float(tokens[5]))
								Y_IMAGE.append(float(tokens[6]))
							
						if 1159<=float(tokens[5])<=1168:                
							if 976<=float(tokens[6])<=986:
								NUMBER.append(float(tokens[0]))
								FLUX_APER.append(float(tokens[1]))
								FLUXERR_APER.append(float(tokens[2]))
								MAG_APER.append(float(tokens[3]))
								MAGERR_APER.append(float(tokens[4]))
								X_IMAGE.append(float(tokens[5]))
								Y_IMAGE.append(float(tokens[6]))
							
				   
				   
						#Comparison Star 1
						if 731<=float(tokens[5])<=740:
							if 1044<=float(tokens[6])<=1054:
								NUMBER_C1.append(float(tokens[0]))
								FLUX_APER_C1.append(float(tokens[1]))
								FLUXERR_APER_C1.append(float(tokens[2]))
								MAG_APER_C1.append(float(tokens[3]))
								MAGERR_APER_C1.append(float(tokens[4]))
								X_IMAGE_C1.append(float(tokens[5]))
								Y_IMAGE_C1.append(float(tokens[6]))
							
						if 1280<=float(tokens[5])<=1290:
							if 967<=float(tokens[6])<=977:
								NUMBER_C1.append(float(tokens[0]))
								FLUX_APER_C1.append(float(tokens[1]))
								FLUXERR_APER_C1.append(float(tokens[2]))
								MAG_APER_C1.append(float(tokens[3]))
								MAGERR_APER_C1.append(float(tokens[4]))
								X_IMAGE_C1.append(float(tokens[5]))
								Y_IMAGE_C1.append(float(tokens[6]))
							
					
						#Comparison Star 2
							
						if 900<=float(tokens[5])<=910:
							if 1225<=float(tokens[6])<=1232:
								NUMBER_C2.append(float(tokens[0]))
								FLUX_APER_C2.append(float(tokens[1]))
								FLUXERR_APER_C2.append(float(tokens[2]))
								MAG_APER_C2.append(float(tokens[3]))
								MAGERR_APER_C2.append(float(tokens[4]))
								X_IMAGE_C2.append(float(tokens[5]))
								Y_IMAGE_C2.append(float(tokens[6]))             
							
						if 1112<=float(tokens[5])<=1125:
							if 788<=float(tokens[6])<=798:
								NUMBER_C2.append(float(tokens[0]))
								FLUX_APER_C2.append(float(tokens[1]))
								FLUXERR_APER_C2.append(float(tokens[2]))
								MAG_APER_C2.append(float(tokens[3]))
								MAGERR_APER_C2.append(float(tokens[4]))
								X_IMAGE_C2.append(float(tokens[5]))
								Y_IMAGE_C2.append(float(tokens[6]))
							

						#Comparison Star 3
							
						if 825<=float(tokens[5])<=831:
							if 926<=float(tokens[6])<=935:
								NUMBER_C3.append(float(tokens[0]))
								FLUX_APER_C3.append(float(tokens[1]))
								FLUXERR_APER_C3.append(float(tokens[2]))
								MAG_APER_C3.append(float(tokens[3]))
								MAGERR_APER_C3.append(float(tokens[4]))
								X_IMAGE_C3.append(float(tokens[5]))
								Y_IMAGE_C3.append(float(tokens[6]))            
					
						if 1188<=float(tokens[5])<=1198:
							if 1084<=float(tokens[6])<=1094:
								NUMBER_C3.append(float(tokens[0]))
								FLUX_APER_C3.append(float(tokens[1]))
								FLUXERR_APER_C3.append(float(tokens[2]))
								MAG_APER_C3.append(float(tokens[3]))
								MAGERR_APER_C3.append(float(tokens[4]))
								X_IMAGE_C3.append(float(tokens[5]))
								Y_IMAGE_C3.append(float(tokens[6]))
							

							
			filein.close()      
							
	"""
	attempt at adding other stars
						#Comparison Star 4
							
						if 547<=float(tokens[5])<=560:
							if 1220<=float(tokens[6])<=1234:
								NUMBER_C4.append(float(tokens[0]))
								FLUX_APER_C4.append(float(tokens[1]))
								FLUXERR_APER_C4.append(float(tokens[2]))
								MAG_APER_C4.append(float(tokens[3]))
								MAGERR_APER_C4.append(float(tokens[4]))
								X_IMAGE_C4.append(float(tokens[5]))
								Y_IMAGE_C4.append(float(tokens[6]))   
							
						if 1455<=float(tokens[5])<=1465:
							if 790<=float(tokens[6])<=800:
								NUMBER_C4.append(float(tokens[0]))
								FLUX_APER_C4.append(float(tokens[1]))
								FLUXERR_APER_C4.append(float(tokens[2]))
								MAG_APER_C4.append(float(tokens[3]))
								MAGERR_APER_C4.append(float(tokens[4]))
								X_IMAGE_C4.append(float(tokens[5]))
								Y_IMAGE_C4.append(float(tokens[6]))  """

							
	#print(FLUX_APER)
	#print(FLUXERR_APER)
	#print(MAG_APER)
	#print(MAGERR_APER)
	#print(NUMBER)
	#print(X_IMAGE)
	#print(Y_IMAGE)
	#print(len(X_IMAGE_C1))
	#print(len(NUMBER_C1))
	#print(len(X_IMAGE_C2))
	#print(len(NUMBER_C2))
	#print(len(X_IMAGE_C3))


	#change lists to arrays
	FLUX_APER=np.array(FLUX_APER)
	FLUXERR_APER=np.array(FLUXERR_APER)
	MAG_APER=np.array(MAG_APER)
	MAGERR_APER=np.array(MAGERR_APER)

	FLUX_APER_C1=np.array(FLUX_APER_C1)
	FLUXERR_APER_C1=np.array(FLUXERR_APER_C1)
	MAG_APER_C1=np.array(MAG_APER_C1)
	MAGERR_APER_C1=np.array(MAGERR_APER_C1)

	FLUX_APER_C2=np.array(FLUX_APER_C2)
	FLUXERR_APER_C2=np.array(FLUXERR_APER_C2)
	MAG_APER_C2=np.array(MAG_APER_C2)
	MAGERR_APER_C2=np.array(MAGERR_APER_C2)

	FLUX_APER_C3=np.array(FLUX_APER_C3)
	FLUXERR_APER_C3=np.array(FLUXERR_APER_C3)
	MAG_APER_C3=np.array(MAG_APER_C3)
	MAGERR_APER_C3=np.array(MAGERR_APER_C3)

	#read in times

	filein = open(tfile, "r")     
	t = []

	for line in filein.readlines():                      #iterates through the lines in the file
				
			#tokens = line.split()                   # breaks line into tokens seperated by ,
			#print(tokens)
		t.append(float(line))
	
	print('read (past) in files')
	
	t=np.array(t)
	
	dataf=[t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3]
	datam=[t,MAG_APER,MAGERR_APER,MAG_APER_C1,MAGERR_APER_C1,MAG_APER_C2,MAGERR_APER_C2,MAG_APER_C3,MAGERR_APER_C3]
	
	return dataf,datam

def magplot(dataf,datam):
	# Plots 
	t,MAG_APER,MAGERR_APER,MAG_APER_C1,MAGERR_APER_C1,MAG_APER_C2,MAGERR_APER_C2,MAG_APER_C3,MAGERR_APER_C3=datam
	t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3=dataf
	
	# plotting magnitudes and flux of wasp
	plt.rcParams["font.family"] = "serif"
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	
	plt.scatter(t,MAG_APER)
	plt.xlabel("time since transit [days]")
	plt.ylabel("Magnitude")
	plt.title('WASP Observed Magnitude')
	plt.show()

	plt.plot(t,FLUX_APER)
	plt.xlabel("time since transit [days]")
	plt.ylabel("Flux")
	plt.title("Wasp Observed Flux")
	#plt.axvline(-0.0274,color='r')  #telescope flip
	#plt.axvline(-0.0297,color='y')  #start transit
	#plt.axvline(0.0377,color='y')#end transit
	#plt.axvline(-0.0297,color='k')  #45 sec exposure
	#plt.axvline(0.0,color='k')  #30 sec exposure
	plt.show()
	
	# plotting fluxes and comparison stars of comparison stars next to WASP
	plt.plot(t,MAG_APER,label='WASP-2')
	plt.plot(t,MAG_APER_C1,label='C1')
	plt.plot(t,MAG_APER_C2,label='C2')
	plt.plot(t,MAG_APER_C3,label='C3')

	#plt.axvline(-0.0274,color='r') #telescope flip
	#plt.axvline(-0.0297,color='y') #start transit
	#plt.axvline(0.0377,color='y') #end transit
	plt.xlabel("time since transit [days]")
	plt.ylabel("Magnitude")
	plt.legend()
	plt.title("Field Star Magnitudes")
	plt.show()

def oFluxPlot(data):
	
	# oFluxPlot(data) plots observed flux during observation (counts per second from scaled images) 
	t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3=data
	plt.errorbar(t,FLUX_APER,yerr=FLUXERR_APER,ecolor='k',label='Wasp2')
	plt.errorbar(t,FLUX_APER_C1,yerr=FLUXERR_APER_C1,ecolor='orange', label='$C_1$')
	plt.errorbar(t,FLUX_APER_C2,yerr=FLUXERR_APER_C2,ecolor='g', label='$C_2$')
	plt.errorbar(t,FLUX_APER_C3,yerr=FLUXERR_APER_C3,ecolor='r', label='$C_3$')

	plt.xlabel("time since transit [days]")
	plt.ylabel("Flux (counts s$^{-1}$)")
	plt.axvline(-0.0377,color='y')
	plt.axvline(0.0377,color='y')
	plt.title("Observed flux of transit observation")
	plt.legend()
	plt.show()

def const_flux_plot(dataf,datam):
	# const_flux_plot(dataf,datam) see if the comparison stars flux and magnitudes subtracted really give straight lines
	t,MAG_APER,MAGERR_APER,MAG_APER_C1,MAGERR_APER_C1,MAG_APER_C2,MAGERR_APER_C2,MAG_APER_C3,MAGERR_APER_C3=datam
	t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3=dataf
	
	MAG_12 = MAG_APER_C1-MAG_APER_C2
	plt.plot(t,MAG_12)
	MAG_13 = MAG_APER_C1-MAG_APER_C3
	plt.plot(t,MAG_13)
	MAG_23 = MAG_APER_C2-MAG_APER_C3
	plt.plot(t,MAG_23)
	MAG_Extinc = MAG_APER-MAG_APER_C1
	plt.plot(t, MAG_Extinc)
	plt.title("Mag subtract")
	plt.show()


	FLUX_12 = FLUX_APER_C1/FLUX_APER_C2
	err_12 = div_err(FLUX_APER_C1,FLUX_APER_C2,FLUXERR_APER_C1,FLUXERR_APER_C2)
	FLUX_13 = FLUX_APER_C1/FLUX_APER_C3
	FLUX_23 = FLUX_APER_C2/FLUX_APER_C3
	FLUX_21 = FLUX_APER_C2/FLUX_APER_C1
	FLUX_31 = FLUX_APER_C3/FLUX_APER_C1
	FLUX_32 = FLUX_APER_C3/FLUX_APER_C2

	plt.errorbar(t,FLUX_12, yerr=err_12, ecolor='k', label='$C_1/C_2 w err$')
	plt.plot(t,FLUX_13, label='$C_1/C_3$')
	plt.plot(t,FLUX_23, label='$C_2/C_3$')
	plt.plot(t,FLUX_23, label='$C_2/C_1$')
	plt.plot(t,FLUX_31, label='$C_3/C_1$')
	plt.plot(t,FLUX_32, label='$C_3/C_2$')
	plt.xlabel("time since transit [days]")
	plt.ylabel("Flux ratio")
	plt.title("Flux ratios of comparison stars")
	plt.legend()
	plt.show()

def extinction_test(data):
	# extinction_test(data) TEST to see how wasp flux reacts to division by comp stars
	t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3=data
	
	FLUX_Extinc = FLUX_APER/FLUX_APER_C1
	FLUX_Extinc2 = FLUX_APER/FLUX_APER_C2
	FLUX_Extinc3 = FLUX_APER/FLUX_APER_C3

	plt.plot(t, FLUX_Extinc, label='$W/C_1$')
	plt.plot(t, FLUX_Extinc2, label='$W/C_2$')
	plt.plot(t, FLUX_Extinc3, label='$W/C_3$')
	plt.title("Flux ratio of WASP/Comparison Stars")
	plt.xlabel("time since transit [days]")
	plt.ylabel("Flux ratio")
	plt.legend()
	plt.show()  


def AtEx(data):
	# determine the mean of comparison star fluxes and divide the flux by mean to get the atmospheric fluctuations from the mean
	# find the weighted average flux fluctuation from all comparison stars
	# WASP flux adjusted for extinction
	# Determine error on WASP as combo of observed wasp errors and reduced flux fluct errors
	# find average flux of WASP2 not during transit to use for normalisation
	# normalise WASP flux and compute errors
	
	t,FLUX_APER,FLUXERR_APER,FLUX_APER_C1,FLUXERR_APER_C1,FLUX_APER_C2,FLUXERR_APER_C2,FLUX_APER_C3,FLUXERR_APER_C3=data
	
	C1_mean_FLUX = np.mean(FLUX_APER_C1)
	#print(np.std(FLUX_APER_C1)/np.sqrt(len(FLUX_APER_C1)))	# error on the mean is 1.6 << FLUXERR so ignored
	C1_Fluct_FLUX = FLUX_APER_C1 / C1_mean_FLUX 			# to do: calc error on the mean, standard deviation
	c1err= div_err(FLUX_APER_C1,C1_mean_FLUX,FLUXERR_APER_C1,0)

	C2_mean_FLUX = np.mean(FLUX_APER_C2)
	C2_Fluct_FLUX = FLUX_APER_C2 / C2_mean_FLUX 
	c2err= div_err(FLUX_APER_C2,C2_mean_FLUX,FLUXERR_APER_C2,0)

	C3_mean_FLUX = np.mean(FLUX_APER_C3)
	C3_Fluct_FLUX = FLUX_APER_C3 / C3_mean_FLUX 
	c3err= div_err(FLUX_APER_C3,C3_mean_FLUX,FLUXERR_APER_C3,0)
	
	# find the weighted mean flux fluctuation (and error) from all comparison stars
	AVG_FLUX_FLUCT,avgerr = w_mean(C1_Fluct_FLUX,C2_Fluct_FLUX,C3_Fluct_FLUX,c1err,c2err,c3err)
	
	# WASP flux adjusted for extinction
	Flux_WASP=FLUX_APER/AVG_FLUX_FLUCT

	# Determine error on WASP as combo of observed wasp errors and reduced flux fluct errors
	fwerr=div_err(FLUX_APER,AVG_FLUX_FLUCT,FLUXERR_APER,avgerr)
	
	# find average flux of WASP2 not during transit to use for normalisation
	# before transit values 0:9, 140:178
	
	index = np.arange(10,139,step=1)   #transit images 10-139 for maybe data, 6-129 for sure data
	WASP_outsidetransit = np.mean(np.delete(Flux_WASP,index))	# calculate std for mean error
	#print(np.std(np.delete(Flux_WASP,index))/np.sqrt(len(np.delete(Flux_WASP,index)))) # error on the mean is 1.6 << fwerr so ignored

	# normalise WASP flux and compute errors
	Flux_WASP_normal = Flux_WASP/WASP_outsidetransit
	fwnerr = (fwerr/Flux_WASP)*Flux_WASP_normal

	'''
	plt.rcParams["font.family"] = "serif"
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	
	plt.plot(t,FLUX_APER_C1,label='$C_1$',color='red')
	plt.plot(t,FLUX_APER_C2, label='$C_2$',color='blue')
	plt.plot(t,FLUX_APER_C3, label='$C_3$',color='darkgrey')
	plt.hlines(C1_mean_FLUX,min(t),max(t),color='darkred')
	plt.hlines(C2_mean_FLUX,min(t),max(t), color='darkblue')
	plt.hlines(C3_mean_FLUX,min(t),max(t),color='k')
	#plt.plot(t, AVG_FLUX_FLUCT)
	plt.title("Observed and Mean Flux of Comparison Field Stars")
	plt.xlabel("Time since mid-transit [days]")
	plt.ylabel("Flux (counts s$^{-1}$)")
	plt.legend()
	plt.show()
	

	plt.plot(t,C1_Fluct_FLUX, label='$C_1$')
	plt.plot(t,C2_Fluct_FLUX, label='$C_2$')
	plt.plot(t,C3_Fluct_FLUX, label='$C_3$')
	plt.plot(t, AVG_FLUX_FLUCT,color='k',label='Mean')
	plt.title("Ratio of flux extinction from comparison field stars")
	plt.xlabel("time since mid-transit [days]")
	plt.ylabel("flux ratio: f$_{obs}/f_{mean}$")
	plt.legend()
	plt.show()
	'''

	#Comparison Stars corrected for extinction, should give straight lines
	Flux_C1=FLUX_APER_C1/AVG_FLUX_FLUCT
	Flux_C2=FLUX_APER_C2/AVG_FLUX_FLUCT
	Flux_C3=FLUX_APER_C3/AVG_FLUX_FLUCT


	plt.rcParams["font.family"] = "serif"
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	
	plt.plot(t,Flux_C1,label='$C_1$ corr',lw=1)
	plt.plot(t,Flux_C2,label='$C_2$ corr',lw=1)
	plt.plot(t,Flux_C3,label='$C_3$ corr',lw=1)
	#plt.plot(t,FLUX_APER_C1, label='$C_1$ obs',color='red')
	#plt.plot(t,FLUX_APER_C2, label='$C_2$ obs',color='blue')
	#plt.plot(t,FLUX_APER_C3, label='$C_3$ obs',color='k')
	#plt.axvline(-0.0377,color='y')  #start transit
	#plt.axvline(0.0377,color='y')#end transit
	plt.title("Comparison Star Atmos.Corrected Fluxes")
	plt.xlabel("Time since mid-transit [days]")
	plt.ylabel("Flux (counts s$^{-1}$)")
	plt.legend()
	plt.show()


	plt.rcParams["font.family"] = "serif"
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	
	#plt.axvline(-0.0274,color='r')  #telescope flip
	#plt.axvline(-0.0377,color='y')  #start transit
	#plt.axvline(0.0377,color='y')#end transit
	#plt.axvline(-0.0297,color='k')  #45 sec exposure
	#plt.axvline(0.0,color='k')  #30 sec exposure
	plt.plot(t,FLUX_APER,label='observed',color='k',lw=1)
	plt.plot(t,Flux_WASP,label='corrected',color='red',lw=1)
	plt.title("Wasp Observed Flux vs Wasp Corrected Flux")
	plt.xlabel("time since transit [days]")
	plt.ylabel("Flux (counts s$^{-1}$)")
	plt.ylim(ymin=2950,ymax=3100)
	plt.legend()
	plt.show()
	
	return Flux_WASP_normal,fwnerr


### Plotting the Model ###

# different models for different fixed parameters
def batmanCurveFit0(xdata,rp0,inc0,e0,w0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w =  w0                       #longitude of periastron (in degrees)
    params.limb_dark = "uniform"        #limb darkening model
    params.u = []      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFitee0(xdata,rp0,inc0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "uniform"        #limb darkening model
    params.u = []      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFit1(xdata,rp0,inc0,e0,w0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w =  w0                       #longitude of periastron (in degrees)
    params.limb_dark = "linear"        #limb darkening model
    params.u = [0.525]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFite1(xdata,rp0,inc0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "linear"        #limb darkening model
    params.u = [0.525]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux    

def batmanCurveFit2(xdata,rp0,inc0,e0,w0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w =  w0                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFite2(xdata,rp0,inc0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux

def batmanCurveFit3(xdata,rp0,inc0,e0,w0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w =  w0                       #longitude of periastron (in degrees)
    params.limb_dark = "nonlinear"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFite3(xdata,rp0,inc0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "nonlinear"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux

def batmanCurveFitAllFree(xdata,rp0,inc0,e0,w0,T0,a0):
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = T0                       #orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = a0                        #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w = w0                        #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFitTmax(xdata,rp0,inc0):	# Tmax , a_mean
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFitTmin(xdata,rp0,inc0,e0,w0):		# Tmin , a_mean
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFitAmax(xdata,rp0,inc0,e0,w0):		# Tmean , a_max
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFitAmin(xdata,rp0,inc0,e0,w0):		# Tmean , a_min
    params = batman.TransitParams()
    
    #rp0,inc0,e0,w0=ig
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    

    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def batmanCurveFite0(xdata,rp0,inc0):
    params = batman.TransitParams()
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = 0.069                       #eccentricity
    params.w =  90.                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    
    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux
    
def AllVarParams(t,Wflux,nWerr):
	# a function to plot the best fit curve allowing all parameters to be free with strict boundary conditions
	ig=[0.135,84.7,0.,0.,2.152060692723822,8.308334998763279]	# final 2, period and semi-major axis, calculated from RV data
	
	s_err=scaleError(nWerr,4.5)	# scale errors
	
	lbound=[0.,0.,0.,0.,2.151805365,8.291457176]			# lower error from RV data
	ubound=[1.,90.,1.,180.,2.152316018,8.37161626187]		# upper error from RV data
	fitThis = spicy.curve_fit(f = batmanCurveFitAllFree, xdata = t, ydata = Wflux, sigma=s_err, p0=ig, bounds=(lbound,ubound))
	
	abest,bbest,cbest,dbest,ebest,fbest = fitThis[0]
	
	aerr = np.sqrt(fitThis[1][0][0])  
	berr = np.sqrt(fitThis[1][1][1])
	cerr = np.sqrt(fitThis[1][2][2]) 
	derr = np.sqrt(fitThis[1][3][3])
	eerr = np.sqrt(fitThis[1][4][4]) 
	ferr = np.sqrt(fitThis[1][5][5])
	#flux=batmanCurveFit3(t,abest,bbest,cbest,dbest)								
	#plt.plot(t, flux, label = ld[3])						# model
		
	print('\033[1m')#+'Limb Darkening model: \'%s\''%ld[3]+'\033[0m')
	print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
	print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
	print('Eccentrcity: %s +/- %s'%(cbest,cerr))
	print('Longitude of periapsis (degrees): %s +/- %s'%(dbest,derr))
	print('Orbital Period (days): %s +/- %s'%(ebest,eerr))
	print('Semi-major axis (a) [R_wasp]: %s +/- %s'%(fbest,ferr))
	print('\033[0m')
	
	params=fitThis[0]
	variance=[aerr,berr,cerr,derr,eerr,ferr]
	
	return params,variance
	

def LDModels(t,data,err,ig):
	# plot figure comparing limb darkening models

	#ig=np.array([0.13,0,0.8,0]) # Rp (R_wasp), orbital inclination, eccentricity, longitude of periapsis
	
	s_err=scaleError(err,4.5)
	
	ld = ["uniform", "linear", "quadratic"]#, "nonlinear"]
	funcs=[batmanCurveFit0,batmanCurveFit1,batmanCurveFit2] #,batmanCurveFit3]


	plt.figure()
	plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)

	
	##### uniform fit #####
	for i in range(len(funcs)):
		fitThis=spicy.curve_fit(f = funcs[i], xdata = t, ydata = data, sigma=s_err, p0=ig, bounds=(0,[1.,90.,1.,360.]))
		abest, bbest,cbest,dbest = fitThis[0]
		aerr = np.sqrt(fitThis[1][0][0])  
		berr = np.sqrt(fitThis[1][1][1])
		cerr = np.sqrt(fitThis[1][2][2]) 
		derr = np.sqrt(fitThis[1][3][3])      
		flux=funcs[i](t,abest,bbest,cbest,dbest)								
		plt.plot(t, flux, label = ld[i])						# model
		
		print('\033[1m'+'Limb Darkening model: \'%s\''%ld[i]+'\033[0m')
		print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
		print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
		print('Ellipticity: %s +/- %s'%(cbest,cerr))
		print('Longitude of periapsis (degrees): %s %s'%(dbest,derr))

	#plt.ylim(0.98, 1.01)
	plt.xlabel("time since transit [days]")
	plt.ylabel("relative flux")
	plt.legend()        
	plt.show()
	
def LDeModels(t,data,err,ig):
	# plot figure comparing limb darkening models

	#ig=np.array([0.13,0,0.8,0]) # Rp (R_wasp), orbital inclination, eccentricity, longitude of periapsis
	
	s_err=scaleError(err,4.5)
	
	ld = ["uniform", "linear", "quadratic"] #, "nonlinear"]
	ld_coefficients = [[], [0.525], [0.525, 0.188]]
	funcs=[batmanCurveFitee0,batmanCurveFite1,batmanCurveFite2]#,batmanCurveFite3]
	cs=["red","lawngreen","k"]

	plt.rcParams["font.family"] = "serif"
	plt.figure()
	plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)

	for i in range(len(funcs)):
	##### uniform fit, linear fit, quadratic fit, nonlinear fit #####
		fitThis=spicy.curve_fit(f = funcs[i], xdata = t, ydata = data, sigma=s_err, p0=ig, bounds=(0,[1.,90.]))
		abest, bbest = fitThis[0]
		aerr = np.sqrt(fitThis[1][0][0])  
		berr = np.sqrt(fitThis[1][1][1])
		#cerr = np.sqrt(fitThis[1][2][2])      
		flux=funcs[i](t,abest,bbest)								

		
		print('\033[1m'+'Limb Darkening model: \'%s\''%ld[i]+'\033[0m')
		print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
		print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
		#print('Longitude of periapsis (degrees): %s %s'%(cbest,cerr))
		rcs=chsq(data,flux,s_err,[abest,bbest])
		plt.plot(t, flux, color=cs[i], label = ld[i]+', $\chi_\\nu^2$=%.3f'%rcs,lw=1)

	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	plt.title('WASP-2: Transit light curve during WASP-2b transit')
	plt.xlabel("Time since mid-transit [days]")
	plt.ylabel("Relative Flux")
	plt.legend()#bbox_to_anchor=(1, 1), loc='upper left')        
	plt.show()
	
def MinChiSqModel(x,y,err,ig):
	params_list=[]
	var_list=[]
	cs_list=[]
	rcs_list=[]
	flux_list=[]
	func_list=[batmanCurveFitTmax,batmanCurveFitTmin,batmanCurveFitAmax,batmanCurveFitAmin]
	s_err=scaleError(err,4.5)					# scales errors by factor of 1/4.5
	
	for i in func_list:
		params,var=findParams(x,y,err,ig,i)
		params_list.append(params)
		var_list.append(var)
		flux=i(x,params[0],params[1],params[2],params[3])
	
		chs,rchs=chsq(y,flux,s_err,params)
		cs_list.append(chs)
		rcs_list.append(rchs)
		
	print(cs_list.index(min(cs_list)))
	bf=cs_list.index(min(cs_list))
	
	print(func_list[bf])
	print(cs_list)
	print(cs_list[bf])
	
	abest,bbest,cbest,dbest=params_list[bf]
	aerr,berr,cerr,derr=var_list[bf]
	
	print('\033[1m')#+'Limb Darkening model: \'%s\''%ld[3]+'\033[0m')
	print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
	print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
	print('Eccentrcity: %s +/- %s'%(cbest,cerr))
	print('Longitude of periapsis (degrees): %s +/- %s'%(dbest,derr))
	print('RChi-squ = %s'%rcs_list[bf])
	print('\033[0m')
	
	
	
	
#### one plot ####

def batmanCurveFit(xdata,rp0,inc0,e0,w0):
    params = batman.TransitParams()
    
    params.t0 = 0.                        #time of inferior conjunction
    params.per = 2.148          			#orbital period days
    params.rp = rp0                       #planet radius (in units of stellar radii)
    params.a = 8.296                  #semi-major axis (in units of stellar radii)
    params.inc = inc0                      #orbital inclination (in degrees)
    params.ecc = e0                       #eccentricity
    params.w =  w0                       #longitude of periastron (in degrees)
    params.limb_dark = "quadratic"        #limb darkening model
    params.u = [0.525, 0.188]      #limb darkening coefficients [u1, u2]
    
    m = batman.TransitModel(params, xdata)
    flux = m.light_curve(params)
    return flux

def findParams(t,W_Flux,nWerr,ig,func):
	
	s_err=scaleError(nWerr,4.5)
	fitThis = spicy.curve_fit(f = func, xdata = t, ydata = W_Flux, sigma=s_err, p0=ig, bounds=(0,[1.,90.,1.,360.]))
	abest,bbest,cbest,dbest = fitThis[0]
	
	aerr = np.sqrt(fitThis[1][0][0])  
	berr = np.sqrt(fitThis[1][1][1])
	cerr = np.sqrt(fitThis[1][2][2]) 
	derr = np.sqrt(fitThis[1][3][3])       
	#flux=batmanCurveFit3(t,abest,bbest,cbest,dbest)								
	#plt.plot(t, flux, label = ld[3])						# model
		
	rj,rjerr=rJup(abest,aerr)
	
	print('\033[1m')#+'Limb Darkening model: \'%s\''%ld[3]+'\033[0m')
	print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
	print('Planet Radius (Jupiter radius): %s +/- %s'%(rj,rjerr))
	print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
	print('Eccentrcity: %s +/- %s'%(cbest,cerr))
	print('Longitude of periapsis (degrees): %s +/- %s'%(dbest,derr))
	print('\033[0m')
	
	params=fitThis[0]
	variance=[aerr,berr,cerr,derr]
	
	return params,variance
	
def findParamse0(t,W_Flux,nWerr,ig,func):
	
	s_err=scaleError(nWerr,4.5)
	fitThis = spicy.curve_fit(f = func, xdata = t, ydata = W_Flux, sigma=s_err, p0=ig, bounds=(0,[1.,90.]))
	abest,bbest = fitThis[0]
	
	aerr = np.sqrt(fitThis[1][0][0])  
	berr = np.sqrt(fitThis[1][1][1])
	#cerr = np.sqrt(fitThis[1][2][2])   
	#flux=batmanCurveFit3(t,abest,bbest,cbest,dbest)								
	#plt.plot(t, flux, label = ld[3])						# model
	
	rj,rjerr=rJup(abest,aerr)
		
	print('\033[1m')#+'Limb Darkening model: \'%s\''%ld[3]+'\033[0m')
	print('Planet Radius (stellar radius): %s +/- %s'%(abest,aerr))
	print('Planet Radius (Jupiter radius): %s +/- %s'%(rj,rjerr))
	print('Oribital Inclination (degrees): %s +/- %s'%(bbest,berr))
	#print('Longitude of periapsis (degrees): %s +/- %s'%(cbest,cerr))
	print('\033[0m')
	
	params=fitThis[0]
	variance=[aerr,berr]
	
	return params,variance


def plotLC(t,data,err,params):
	# plot light curve for 4 best fit parameters for Rp, incl, ecc, LoP
	abest,bbest,cbest,dbest=params
	flux= batmanCurveFit(t,abest,bbest,cbest,dbest)
	s_err=scaleError(err,4.5)					# scales errors by factor of 1/4.5
	
	plt.rcParams["font.family"] = "serif"
	plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)
	plt.plot(t, flux, color='darkorange',zorder=2,lw=1)
	
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	plt.title('WASP-2: Transit light curve during WASP-2b transit')
	plt.xlabel("Time since transit [days]")
	plt.ylabel("Relative Flux")
	plt.show()
	
def plotLCe0(t,data,err,params):
	# plot light curve for 4 best fit parameters for Rp, incl, ecc, LoP
	abest,bbest=params
	flux= batmanCurveFite0(t,abest,bbest)
	s_err=scaleError(err,4.5)					# scales errors by factor of 1/4.5
	
	chsq(data,flux,s_err,params)
	
	model=str(input("Fit model curve? yes or no: "))
	plt.rcParams["font.family"] = "serif"
	if model=='no':
		plt.scatter(t, data ,color='blue',zorder=1,marker='.',s=5)
	else:
		plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)
		plt.plot(t, flux, color='darkorange',zorder=2,lw=1)
	
	#scale_err(data,flux,err) #tests by what factor to scale down errors such that 2/3 cross expected value
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	#plt.title('WASP-2: Transit light curve during WASP-2b transit')
	plt.title('WASP-2: Processed flux data during WASP-2b transit')
	plt.xlabel("Time since mid-transit [days]")
	plt.ylabel("Relative Flux")
	plt.show()

def plotLCall(t,data,err,params):
	# plot light curve for 6 best fit parameters for Rp, incl, ecc, LoP, T and a
	abest,bbest,cbest,dbest,ebest,fbest=params
	flux= batmanCurveFitAllFree(t,abest,bbest,cbest,dbest,ebest,fbest)
	s_err=scaleError(err,4.5)					# scales errors by factor of 1/4.5
	
	plt.rcParams["font.family"] = "serif"
	plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.1)
	plt.plot(t, flux, color='darkorange',zorder=2)
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	plt.title('WASP-2: Transit Light Curve of Wasp-2b $_{all~params~free}$')
	plt.xlabel("Time since transit [days]")
	plt.ylabel("Relative Flux")
	plt.show()

### Plotting the LC + residual plots ###

def residuals(t, nWflux, nWerr, params):
	'''
	t 		- time data
	Wflux 	- normalised flux
	s_err	- scaled errors
	params 	- best fit params
	'''
	abest,bbest,cbest,dbest=params
	flux= batmanCurveFit(t,abest,bbest,cbest,dbest)
	s_err=scaleError(nWerr,4.5)					# scales errors by factor of 1/4.5
	
	chsq(nWflux,flux,s_err,params)
	
	plt.rcParams["font.family"] = "serif"
	
	fig,(ax1,ax2)=plt.subplots(2,sharex=True)

	ax1.xaxis.set_minor_locator(AutoMinorLocator())
	ax1.yaxis.set_minor_locator(AutoMinorLocator())
	ax1.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	ax1.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	ax2.xaxis.set_minor_locator(AutoMinorLocator())
	ax2.yaxis.set_minor_locator(AutoMinorLocator())
	ax2.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	ax2.tick_params(which='minor', length=4, direction='in', top=True, right=True)

	fig.suptitle('WASP-2: Transit light curve during WASP-2b transit')
	ax1.errorbar(t, nWflux ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)
	ax1.plot(t, flux, color='darkorange',zorder=2,lw=1)
	ax2.axhline(0,ls='-',alpha=1,c='darkorange',zorder=1,lw=1)
	#print((t),min(t),max(t))
	ax2.errorbar(t,(nWflux-flux),yerr=s_err,fmt='.',ms=1.5,ecolor='k',color='blue',elinewidth=0.15)
	#plt.axvline(-0.0377,color='y')
	#plt.axvline(0.0377,color='y')
	#plt.ylim(0.96, 1.02)
	ax2.set_xlabel("Time since mid-transit [days]")
	ax2.set_ylabel("Relative Flux Residuals")
	ax1.set_ylabel("Relative Flux")
	#ax1.legend()        
	plt.show()
	
def residualse0(t, nWflux, nWerr, params):
	'''
	t 		- time data
	Wflux 	- normalised flux
	s_err	- scaled errors
	params 	- best fit params
	'''
	abest,bbest=params
	flux= batmanCurveFite0(t,abest,bbest)
	s_err=scaleError(nWerr,4.5)					# scales errors by factor of 1/4.5
	
	chsq(nWflux,flux,s_err,params)
	
	plt.rcParams["font.family"] = "serif"
	
	fig,(ax1,ax2)=plt.subplots(2,sharex=True)

	ax1.xaxis.set_minor_locator(AutoMinorLocator())
	ax1.yaxis.set_minor_locator(AutoMinorLocator())
	ax1.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	ax1.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	ax2.xaxis.set_minor_locator(AutoMinorLocator())
	ax2.yaxis.set_minor_locator(AutoMinorLocator())
	ax2.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	ax2.tick_params(which='minor', length=4, direction='in', top=True, right=True)

	fig.suptitle('WASP-2: Transit light curve during WASP-2b transit')
	ax1.errorbar(t, nWflux ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.15)
	ax1.plot(t, flux, color='darkorange',zorder=2,lw=1)
	ax2.axhline(0,ls='-',alpha=1,c='darkorange',zorder=1,lw=1)
	#print((t),min(t),max(t))
	ax2.errorbar(t,(nWflux-flux),yerr=s_err,fmt='.',ms=1.5,ecolor='k',color='blue',elinewidth=0.15)
	#plt.axvline(-0.0377,color='y')
	#plt.axvline(0.0377,color='y')
	#plt.ylim(0.96, 1.02)
	ax2.set_xlabel("Time since mid-transit [days]")
	ax2.set_ylabel("Relative Flux Residuals")
	ax1.set_ylabel("Relative Flux")
	#ax1.legend()        
	plt.show()


def modelcomp(t,data,err,m1,m2):
	s_err=scaleError(err,4.5)
	rchs1=chsq(data,m1,s_err,[1,2,3,4])
	rchs2=chsq(data,m2,s_err,[1,2])
	
	plt.rcParams["font.family"] = "serif"
	plt.errorbar(t, data ,yerr=s_err,color='blue',ecolor='k',zorder=1,fmt='.',ms=1.5,elinewidth=0.2)
	plt.plot(t, m1, color='purple',zorder=2,label='$e=0.35(23)$\n $\chi_\\nu^2$=%.3f'%rchs1)
	plt.plot(t, m2, color='darkorange',zorder=3,label='$e=0.07(9)$\n $\chi_\\nu^2$=%.3f'%rchs2)
	
	
	plt.minorticks_on()
	plt.tick_params(which='major', direction='in', length=6, width=2, colors='black', top=True, right=True)
	plt.tick_params(which='minor', length=4, direction='in', top=True, right=True)
	plt.title('WASP-2: Transit Light Curve of WASP-2b')
	plt.xlabel("Time since transit [days]")
	plt.ylabel("Relative Flux")
	plt.legend()
	plt.show()

def main():

	dataf,datam=read()
	t=dataf[0]
	#magplot(dataf,datam)			# plots wasp mag and flux and comp star mag
	#oFluxPlot(dataf) 				# plots observed flux and raw errs of all stars during observation (counts per second from scaled images)
	#const_flux_plot(dataf,datam) 	# see if the comparison stars flux and magnitudes subtracted really give straight lines 
	#extinction_test(dataf) 		# TEST to see how wasp flux reacts to division by comp stars
	nWflux,nWerr=AtEx(dataf)
	
	
	ig=[0.135,84.7,0.,0.] # Rp,incl,eccent,LoP
	ige0=[0.135,84.7] # Rp, incl, fixing e and LoP from RV analysis
	
	LDeModels(t,nWflux,nWerr,ige0)    										# determines quadratic limb darkening fit is best
	#MinChiSqModel(t,nWflux,nWerr,ig)										# investigates uncertainty variation of free parameters by varying fixed parameters
	
	params,var=findParams(t,nWflux,nWerr,ig,batmanCurveFit)					# finds 4 free parameters (incl. eccentricity)
	paramse0,vare0=findParamse0(t,nWflux,nWerr,ige0,batmanCurveFite0)		# finds 3 free parameters (excl. eccentricity)
	m1=batmanCurveFit(t,params[0],params[1],params[2],params[3])			# returns model for free eccentricity
	m2=batmanCurveFite0(t,paramse0[0],paramse0[1])							# returns model for fixed eccentricity
	modelcomp(t,nWflux,nWerr,m1,m2)											# plots m1 and m2 comparison with reduced chi-squared test
	
	residualse0(t,nWflux,nWerr,paramse0)					# plots eccentricity fixed curve and residuals
	
	plotLCe0(t,nWflux,nWerr,paramse0)						# plots eccentricity fixed light curve only
	
	#Plots curve for most all parameters free (old, don't use)
	'''
	#allp,allv=AllVarParams(t,nWflux,nWerr)	
	#plotLCall(t,nWflux,nWerr,allp)
	'''

	
main()
	




# actually obital period = 2.15d   rp = 0.1350559527 R_Wasp  ... wikipedia
# semi major axis 8.240983371 R_Wasp


#(array([5.2187809 , 0.50045221]), array([[ 2.43943077e-01, -2.09411105e-03],
#       [-2.09411105e-03,  3.73990499e-05]])) for e=0.4 and inc=80