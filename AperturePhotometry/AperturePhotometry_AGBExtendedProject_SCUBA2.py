from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import matplotlib.gridspec as gridspec
from scipy.integrate import simps
from numpy import trapz
import matplotlib.ticker as mtick
import photutils as photutils
import csv
from astropy.table import Table, Column,  vstack
from astropy.io import ascii

#Object fits file input
def observation_file(input_filename, sq_pix_size, wavelength):
	
	fitsFile = fits.open(input_filename)

	data = fitsFile[0].data[0] #use 0th frame (primary image/frame) of the fits file for SCUBA2
	data[np.isnan(data)] = 0
	data = data*sq_pix_size #Converting data from mJy/arcsec^2 to mJy 	

	return data

#Error fits file input
def observation_file_unc(input_filename, sq_pix_size, wavelength):

	fitsFile_error = fits.open(input_filename)

	data_unc = fitsFile_error['VARIANCE'].data[0] #BUNIT = '(mJy/arcsec^2)**2' - Units of the Variance array.	
	data_unc = (np.sqrt(data_unc))*sq_pix_size #Convert from variance to std.dev with mJy units
	
	return data_unc


#Aperture Photometry
def aperture_photometry_data(star, input_filename, r_stellar, r_inner, r_outer, central_pixel_x, central_pixel_y, wavelength, data, data_unc):

	#Openeing Scuba2 data and error file
	data = observation_file(input_filename, sq_pix_size, wavelength)
	data_unc = observation_file_unc(input_filename, sq_pix_size, wavelength)

	#Getting stellar central postion from fits file header
	fitsFile = fits.open(input_filename)	
	cx = central_pixel_x 
	#print ('cx =', cx)
	cy = central_pixel_y
	#print ('cy =', cy)
	central_positions = (cx-1, cy-1) 

	#Aperture Photometry
	stellar_aperture = photutils.CircularAperture(central_positions, r_stellar) #Stellar Aperture
	sky_annlus = photutils.CircularAnnulus(central_positions, r_inner, r_outer) #Sky Annulus

	error = data_unc
	apertures = [stellar_aperture, sky_annlus]
	aperture_photometry = photutils.aperture_photometry(data, apertures, error=error) #Carrying out Aperture Photometry

	#Adding Object name Aperture Photometry Output Table
	aperture_photometry['Object'] = star

	#x and y centers of the annuli are already automatically in the aperture photometry tables found using Aperture Photometry tool using the measured 3sigma radii. 

	#Calculating Final Mean Stellar Flux and Unc. 
	bkg_mean = aperture_photometry['aperture_sum_1'] / sky_annlus.area() #Caclulating BG flux
	bkg_sum = bkg_mean * stellar_aperture.area()

	stellar_flux =  aperture_photometry['aperture_sum_0'] - bkg_sum #Final Mean Flux (Background reduced)
	aperture_photometry['final_stellar_flux_sum'] = stellar_flux 

	aperture_photometry_unc = np.sqrt(

					((aperture_photometry['aperture_sum_err_0'])**2) + 
					
						(
							((aperture_photometry['aperture_sum_err_1'])**2)
							*(stellar_aperture.area()/sky_annlus.area())
						)
					)


	aperture_photometry['final_stellar_flux_uncertainty'] = aperture_photometry_unc
	
	#print(aperture_photometry)

	return aperture_photometry	




if __name__=="__main__":

	wavelength = '850'
	sq_pix_size = 16 # 450mic=4 | 850mic=16 cause unit=mJy/arcsec^2
	source_list_filename = 'Source_List_' + wavelength + '.csv'
	source_list = Table.read(source_list_filename, format='ascii.csv')
	output_filename = 'AperturePhotometry_' + wavelength + '.dat'

	#exit()

	for source in source_list: 

		star = str(source['Source'])
		central_pixel_x = source['Center_X']
		central_pixel_y = source['Center_Y']
		r_stellar = source['Rstar'] #Object Annulus. #Pixel units 
		r_inner = source['Rin'] #Inner Sky Annulus #Pixel Units 
		r_outer = source['Rout'] #Outer Sky Annulus #Pixel Units 

		input_filename = star + '_' + wavelength + '.fits'

		data = observation_file(input_filename, sq_pix_size, wavelength)
		data_unc = observation_file_unc(input_filename, sq_pix_size, wavelength)
		
		aperture_photometry_output = Table(aperture_photometry_data(star, input_filename, r_stellar, r_inner, r_outer, central_pixel_x, central_pixel_y, wavelength, data, data_unc)) #Table() converts the aperture_photometry_data output into an astropy table
		print(aperture_photometry_output) 

		for row in aperture_photometry_output:

			print row['Object'], row['xcenter'], row['ycenter'], row['final_stellar_flux_sum'], row['final_stellar_flux_uncertainty']
		
			f = open(output_filename, 'a')

			outstring = str(row['Object']) + ","  + str(row['xcenter']) + "," +  str(row['ycenter']) + "," + str(row['final_stellar_flux_sum']) + "," + str(row['final_stellar_flux_uncertainty']) + "\t" + "\n"  
			f.write(outstring)	
			f.close()
