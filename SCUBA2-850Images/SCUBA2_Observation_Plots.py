import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import aplpy
from astropy import units as u
from astropy.table import Table

font = {'family' : 'normal',
        'size'   : 18,
	'style'  : 'oblique',
	'weight' : 'medium'}

plt.rc('font', **font)

"""
###################################################################################################################################
    ########### Radial Profile Generation ################################
###################################################################################################################################

- Script to plot the image of the source at a chosen SCUBA-2 wavelength.

- The PSF shape and three sigma brightness extension is overplotted on the image.

- While the given script is for SCUBA-2 850micron it can be adjusted for 450micron as well.  

** This script can be used for both 450 micron obs and 850 micron obs by changing the wavelength varaible at the top of the script. We have only generated the 850micron obs images for this paper. **


###### Input files required ####################
1) .csv table file containing source central X and Y positions and three sigma radius - Source_Information_SCUBA2_850.csv
2) .fits 850micron SCUBA-2 observation of chosen source. 


######## Output ####################

- Image of the source with overplotted PSF & three sigma brightness extension. Axes provided as RA and Dec in degrees as well as the surface brightness colour scale in mJy/arcsec^2. 




"""



wavelength = 850 #Change to match wavelength. Can be adjusted to generate 450 observations as well. The cetral position, Vmin and Vmax values will have to be adjusted accordingly

stretch_type_dict = {'850':'arcsinh', '450':'arcsinh'} 
stretch_type = stretch_type_dict[str(wavelength)]

cmap_type_dict =  {'850':'magma', '450':'magma'}
cmap_type = cmap_type_dict [str(wavelength)]

beam_fwhm_dict = {'850': 13/3600., '450': 7.9/3600.} #SCUBA2 Beam FWHM. #Units=degrees #850micron main beam =13"=0.00361degrees #450micron main beam =7.9"=0.00219444degrees
beam_fwhm = beam_fwhm_dict[str(wavelength)] 

fig_size_dict = {'850': 90/3600., '450': 70/3600.}
fig_size = fig_size_dict[str(wavelength)]

Source_Data = Table.read('Source_Information_SCUBA2_'+str(wavelength)+'.csv', format='ascii.csv')

for Source in Source_Data: 


	star = Source['Source']
	x_center = Source['Center_X'] #Units in pixels 
	y_center = Source['Center_Y']  #Units in pixels 
	three_sigmaRadius = Source['ThreeSigma_Radius(arcsec)'] #Units =arcsec
	vmax_val = Source['Vmax_use'] #Vmax value estimated by checking the value of the brightest pixel and then reducing a small value from it to allow for saturation. 

	#Creating Plot using Aplpy
	fig= aplpy.FITSFigure(star+'_'+str(wavelength)+'.fits') #Load fits file
	fig.show_colorscale(vmin=0.001, vmax=vmax_val, cmap=cmap_type, stretch=stretch_type) #vmin=0.0001 for log scale - vmin can not be zero for log scale so we chose an arbitrary number close to zero| vmax=vmax_val can be set if we want to.  

	#Cropping plot to the size we want
	x_world, y_world = fig.pixel2world(x_center, y_center) #Converting pixels into degrees for plotting 
	fig.recenter(x_world, y_world, 90/3600.) #Units in degrees || (x, y, size_in_degrees) ||90'' = 0.025degrees

	#Adding Colourbar
	fig.add_colorbar()
	fig.colorbar.set_width(0.3)
	fig.colorbar.set_location('right')
	fig.colorbar.set_axis_label_text('Surface Brightness (mJy/arcsec$^{2}$)')
	fig.colorbar.set_axis_label_pad(20)
	fig.colorbar.set_axis_label_rotation(270)


	beam_radius = beam_fwhm/2 #rad = FWHM/2 (diam=fwhm) #units=degrees
	three_sigmaRadius_degrees = three_sigmaRadius / 3600. #Converting from arcseconds to degrees.  

	#Over-plotting beam and 3 sigma extensions as a circles 
	fig.show_circles(x_world, y_world, beam_radius, edgecolor='darkslategrey', dashes='--', linewidth=2) #Prim.Beam (xcenter, ycenter, radius)
	fig.show_circles(x_world, y_world, three_sigmaRadius_degrees, edgecolor='red', dashes='--', linewidth=2) #3sig.Radius (xcenter, ycenter, radius)

	fig.save(star + '_' + 'ObsImage' + '_' + str(wavelength) + '.png', dpi=300)
	#plt.show()







