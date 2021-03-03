# TGP-Exoplanets
Senior Honours Telescope Group Project - March 2021

This repository contains all of the code we have written to aid us with our Senior Honours Telescope Group Project at The University of Edinburgh. 

The files have been written by Euan Newlands, Gemma Robertson, Lucie Scharre, Douglas Tait, Victoria Thompson and Gregory Wilcox.


MultiParam.py
This file was part of LS and EN's work on the light curve fitting process.
Description: This code reads in time-scaled catalogs from WASP-2 (WASAP-2b transit) observation. 
WASP-2 and 3 comparison field stars are used read in. Atmospheric effects are corrected for by using the field stars
where the final corrected WASP-2 flux is plotted and a light curve is fitted using the batman package.
The best fit light curve is determined and the model returns best fit parameters, used to determine the properties of WASP-2b


dustextinction.py
This file was created by GR and VT for the purpose of  


dustextinction.ipynb
This notebook was created by DT for the purpose of producing graphs for the dust extinction in each filter.


RV_Fitting.ipynb
This notebook is an old method used by DT and GW for the purpose of radial velocity fitting. 
It is now out of use as we changed our model however it has been left in for completeness. 


rv.py
This file was created by DT for the purpose of radial velocity fitting. It built on work by Adam Dempsey who had written the core of the code. (https://adamdempsey90.github.io/python/radial_velocities/radial_velocities.html)
This file is the brains of the program. It creates a star object, with its mass and RV data.
The file then contains all of the functions necessary to calculate and return planetary parameters and create relevant graphs.


rv_fit.py
This file was created by DT for the purpose of radial velocity fitting. It built on work by Adam Dempsey who had written the core of the code. (https://adamdempsey90.github.io/python/radial_velocities/radial_velocities.html)
This file reads in the data and performs the radial velocity fit using scipy.optimize.cure_fit().
It then prints out a list of the fitted and calculated planetary properties.
It also produces graphs of the curve fit, photophase fit and residuals.

