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







