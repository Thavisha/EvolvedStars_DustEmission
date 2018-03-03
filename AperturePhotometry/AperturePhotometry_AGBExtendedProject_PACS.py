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

"""
###################################################################################################################################
    ########### Aperture Photometry on SCUBA-2 and PACS observations using Photutils ################################
###################################################################################################################################

- Carry out mean sky reduced aperture photometry on SCUBA-2 and PACS reduced .fits files in order to get the total source flux.

- The central pixel and sky aperture radii were initally measured using CalTech Aperture Photometry Tool (APT) and are listed in table form.
- The source aperture raidus is set to the extension at 3sigma brightness level (converted to pixel units) measured using the radial profiles and listed in the paper. 

- Two scripts, one for each instrument. Each script can be run for either of the two wavelengths. The wavelength must be specified in the script.  

- Output - .dat table containing i)source, ii)central pixel position X, iii)central pixel position Y, iv)source flux, v)flux unc. for the selected wavelength




###### Input files required ####################

#### .csv Tables containing central pixel and aperture radii information ####
1)Source_List_70.csv
2)Source_List_160.csv
3)Source_List_450.csv
4)Source_List_850.csv 


#### .fits Files ####
1) .fits files of the sources listed in "Source_List_.." table for the chosen wavelength 



"""

#Object fits file input
def observation_file(input_filename, sq_pix_size, wavelength):
	
	fitsFile = fits.open(input_filename)

	data = fitsFile['image'].data
	data = data*sq_pix_size #Converting data from mJy/pixel to mJy 

	return data

#Error fits file input
def observation_file_unc(input_filename, sq_pix_size, wavelength):

	fitsFile_error = fits.open(input_filename)

	data_unc = fitsFile_error['stDev'].data 
	data_unc = data_unc * sq_pix_size 
	
	return data_unc


#Aperture Photometry
def aperture_photometry_data(star, input_filename, r_stellar, r_inner, r_outer, central_pixel_x, central_pixel_y, wavelength, data, data_unc):

	#Openeing data and error file
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

	wavelength = '70'
	sq_pix_size = 1 #70mic=1 & 160mic=1 cause PACS data are given in Jy/pixel | 450mic=4 | 850mic=16 cause unit=mJy/arcsec^2
	source_list_filename = 'Source_List_' + wavelength + '.csv'
	source_list = Table.read(source_list_filename, format='ascii.csv')
	output_filename = 'UAnt_AperturePhotometry_' + wavelength + '.csv'

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
