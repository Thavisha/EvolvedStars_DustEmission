###################################################################################################################################
    ########### Dust and Gass Mass Loss ################################
###################################################################################################################################

- Calculate the total dust and gas mass loss, mass loss rates and dust-to-gass ratios.

- Total dust mass derived by integrating over the density profile resulting from SED fitting. 
- projected size converted to physical size using terminal velocity and distance to source using the values presented in the paper. 
- Total gas mass derived using CO gas mass loss rates from DeBeck et al., 2010 and the age of CSE at the radial point. 
- All values derived up to the 3 sigma extension point of 160micron as we need a minimum of two detections at each radial point in order to successfully derive the density at the given radial point, therefore the more extended 70micron radius is unsuitable. But the script will work for any chosen radius.  


###### Input files required ####################

1) 160_3sigma_Extensions_And_Total_Fluxes_And_DeBeckMassLossRates.csv - table containing 3 sigma extensions, DeBeck data and distances, etc needed. 
2) x_interp.csv - file containing x radial points
3) .csv files containing the density profiles derived from SED fitting section for all source listed in file 1). eg: cit6_density_radial_output.csv


##### Output ############

.csv table containing i)source, ii)3 sigma radius in arcsec, iii)3 sigma radius in cm, iv)3 sigma radius in pc, v)DeBeck mass loss rate, vi) v)DeBeck terminal velocity, vi)age of CSE at 3sigma radius, vii)total gas mass, viii)Dust mass using derived gas mass and 1/200 d:g ratio, ix) Total dust mass from integrating the density profile, x)Ratio of the two dust masses, xi) Dust-to-Gas ratio from ix/vii.   



################################################################################################################################
################## Plotting Gas MLR vs. Dust MLR - with best fit line and Dust-to-Gas ratio = 1/200 line ########################
#################################################################################################################################

- Script to plot the CO Gas MLR vs. Dust MLR (dervied from integrating the density profile) with Chi sqaure minimized line of best fit and the canonical dust-to-gas ratio of 1/200 line. 

- Plotting results derived from the script described above. 

###### Input files required ####################
1) GasMLR_vs_DustMLR_withAccepted+BestFitLines.csv - table containing the information required to produce the plot. 


########## Output ##############
- plot of  CO Gas MLR vs. Dust MLR with Chi sqaure minimized line of best fit and the canonical dust-to-gas ratios for C and O rich.



