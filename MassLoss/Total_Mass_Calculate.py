import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from scipy.integrate import simps
from numpy import trapz
import matplotlib.ticker as mtick
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda, blackbody_nu
from astropy.table import Table
from astropy.constants import c #c=speed of light 

"""
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
"""


#Total Mass from the Density Profile Derived from SEDs. #Integration of Density Profile
def total_mass_from_density_profile(star, three_sigma_cm, distance, nskip=1): #nskip=1 - Skip values in the 0th radial point

	density_data_filename = star + '_density_radial_output.csv'
	density_data = Table.read(density_data_filename, format='ascii.csv')

	density_profile = density_data['Density']
	density_profile = density_profile[nskip:]
	density_profile_plusUnc = density_data['Plus_Unc']
	density_profile_plusUnc = density_profile_plusUnc[nskip:]
	density_profile_minusUnc = density_data['Minus_Unc']
	density_profile_minusUnc = density_profile_minusUnc[nskip:]

	x = density_data['x_interp'] #In arcsec
	x = x[nskip:]
	x_cm = (distance * x) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	three_sigma_cm_limit = x_cm < three_sigma_cm
	
	#make a new axis staring from min x interp to max interp in units of 1???

	Integrated_Density = simps((2*np.pi*x_cm[three_sigma_cm_limit]*(10**density_profile[three_sigma_cm_limit])), x_cm[three_sigma_cm_limit]) #Units = g
	Integrated_Density = Integrated_Density /  1.9884754e33 #Convert to Solar Masses
	print ('Integrated_Density =', Integrated_Density) 

	#uncertainty ???

	return Integrated_Density



if __name__=="__main__":

	wavelength = '160' #Get masses for the 160 3 sigma cause you need at least two radial points to get an accurate density. so 70 extension is not enough. at 160 there's both 70 and 160 detections
	wavelength_val = 160
	data_filename = wavelength + '_3sigma_Extensions_And_Total_Fluxes_And_DeBeckMassLossRates.csv'
	input_data = Table.read(data_filename, format='ascii.csv')
	print input_data.columns
	calculated_densitites_output_filename = '000_TotalMass_output_' + wavelength + '.csv'
	#print source_list


	f = open(calculated_densitites_output_filename, 'w')
	f.write("source" + "," + wavelength+"_three_sigma(arcsec)" + "," + wavelength+"_three_sigma(cm)" + "," + wavelength+"_three_sigma(pc)" + "," + "DeBeck_CO_Mass_Loss_Rate (M_sun/yr)" + "," + "DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)" + "," + "3sigma_radius_160_Time(years)(Time=Dist/Vel)unit_sToYr=*3.171e-8" + "," + "DeBeck_Total_CO_Mass_loss_In_3sigma_time_(Msun)" + "," + "Dust_Mass_from_d:g=1:200_ratio(Msun)_(Maybe only the central component??)" + "," + "Intergrated_Dust_Mass_fromDensityProfile(Msun)" + "," + "Integrated_DustMass/DeBeck_DustMass" + "," + "Dust:Gas_Ratio_IntegratedDustMass/DeBeckCOMass_NormalAccepted=1/200_(0.005)" + "\n")
	f.close()

	for source in input_data: 

		star = str(source['Source'])
		distance = source['Distance(pc)'] #In pc
		distance_cm = distance * 3.086E+18 #Convert from pc to cm
		three_sigma_arcsec = source['three_sigma_Extension(arcsec)'] #In arcsec. needs to be converted to cm using the source distances
		three_sigma_cm = (distance * three_sigma_arcsec) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
		three_sigma_pc = three_sigma_cm * 3.24078E-019 
		DeBeck_MassLoss_Rate =  source['DeBeck_CO_Mass_Loss_Rate(M_sun/yr)']
		DeBeck_Terminal_Velocity = source['DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)']

		
		Integrated_Density = total_mass_from_density_profile(star, three_sigma_cm, distance, nskip=1)

		three_sigma_radius_inTimeUnits = (three_sigma_cm/DeBeck_Terminal_Velocity)*3.171E-8 #3sigma_radius_160_Time (years) (Time=Dist/Vel) unit_s to Yr = *3.171e-8. Same as doing (seventy_3sigma_cm / velocity) / 3.154e+7 #/3.154e+7 to convert from s to years

		DeBeck_total_CO_MassLoss_in_3sigma_time = DeBeck_MassLoss_Rate * three_sigma_radius_inTimeUnits

		Dust_Mass_fromDeBeck = DeBeck_total_CO_MassLoss_in_3sigma_time * 0.005 #Accepted dust:gas ratio = 1/200 = 0.005

		Dust_Mass_Ratio = Integrated_Density / Dust_Mass_fromDeBeck 

		Dust_to_Gas_Ratio = Integrated_Density / DeBeck_total_CO_MassLoss_in_3sigma_time


		
		f = open(calculated_densitites_output_filename, 'a')
		outstring = str(star)  + ","  + str(three_sigma_arcsec)  + ","  + str(three_sigma_cm)   + ","  + str(three_sigma_pc) + "," + str(DeBeck_MassLoss_Rate) + "," + str(DeBeck_Terminal_Velocity) + "," + str(three_sigma_radius_inTimeUnits) + "," + str(DeBeck_total_CO_MassLoss_in_3sigma_time) + "," + str(Dust_Mass_fromDeBeck) + "," + str(Integrated_Density) + "," + str(Dust_Mass_Ratio) + "," + str(Dust_to_Gas_Ratio) + "\n"
		f.write(outstring)	
		f.close()


		























