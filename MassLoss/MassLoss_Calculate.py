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


#Total Mass from the Density Profile Derived from SEDs. #Integration of Density Profile
def total_mass_from_density_profile(star, three_sigma_cm, distance):

	density_data_filename = star + '_density_radial_output.csv'
	density_data = Table.read(density_data_filename, format='ascii.csv')
	x = density_data['x_interp'] #In arcsec
	x_cm = (distance * x) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	density_profile = density_data['Density']
	density_profile_plusUnc = density_data['Plus_Unc']
	density_profile_minusUnc = density_data['Minus_Unc']
	three_sigma_cm_limit = x_cm < three_sigma_cm

	Integrated_Density = simps((2*np.pi*x_cm[three_sigma_cm_limit]*(10**density_profile[three_sigma_cm_limit])), x_cm[three_sigma_cm_limit]) #Units = g
	Integrated_Density = Integrated_Density /  1.9884754e33 #Convert to Solar Masses
	print ('Integrated_Density =', Integrated_Density) 


	return Integrated_Density



if __name__=="__main__":

	wavelength = '160' #Get masses for the 160 3 sigma because we need detections at least two radial points to get an accurate density. At 160 there's both 70 and 160 detections
	wavelength_val = 160
	data_filename = wavelength + '_3sigma_Extensions_And_Total_Fluxes_And_DeBeckMassLossRates.csv'
	input_data = Table.read(data_filename, format='ascii.csv')
	print input_data.columns
	calculated_densitites_output_filename = 'TotalMass_output_' + wavelength + '.csv'
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
		three_sigma_pc = three_sigma_cm *3.24078E-019
		total_flux = source['Total_Flux(Jy)']
		total_flux_unc = source['Total_Flux_Unc(Jy)']  
		DeBeck_MassLoss_Rate =  source['DeBeck_CO_Mass_Loss_Rate(M_sun/yr)']
		DeBeck_Terminal_Velocity = source['DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)']

		
		Integrated_Density = total_mass_from_density_profile(star, three_sigma_cm, distance)

		three_sigma_radius_inTimeUnits = (three_sigma_cm/DeBeck_Terminal_Velocity)*3.171E-8 #3sigma_radius_160_Time (years) (Time=Dist/Vel) unit_s to Yr = *3.171e-8

		DeBeck_total_CO_MassLoss_in_3sigma_time = DeBeck_MassLoss_Rate * three_sigma_radius_inTimeUnits

		Dust_Mass_fromDeBeck = DeBeck_total_CO_MassLoss_in_3sigma_time * 0.005 #Accepted dust:gas ratio = 1/200 = 0.005

		Dust_Mass_Ratio = Integrated_Density / Dust_Mass_fromDeBeck 

		Dust_to_Gas_Ratio = Integrated_Density / DeBeck_total_CO_MassLoss_in_3sigma_time


		
		f = open(calculated_densitites_output_filename, 'a')
		outstring = str(star)  + ","  + str(three_sigma_arcsec)  + ","  + str(three_sigma_cm)   + ","  + str(three_sigma_pc) + "," + str(DeBeck_MassLoss_Rate) + "," + str(DeBeck_Terminal_Velocity) + "," + str(three_sigma_radius_inTimeUnits) + "," + str(DeBeck_total_CO_MassLoss_in_3sigma_time) + "," + str(Dust_Mass_fromDeBeck) + "," + str(Integrated_Density) + "," + str(Dust_Mass_Ratio) + "," + str(Dust_to_Gas_Ratio) + "\n"
		f.write(outstring)	
		f.close()


		























