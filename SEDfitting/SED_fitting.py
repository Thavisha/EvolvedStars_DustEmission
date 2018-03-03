import matplotlib.pyplot as plt
import numpy as np
from astropy import constants as const
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from scipy.optimize import curve_fit
from astropy.table import Table
from astropy.constants import c #c=speed of light 
from numpy import loadtxt
from scipy import interpolate
import emcee
import sys
import corner
from scipy.stats import norm
from numpy import where 
from numpy import asarray
from astropy.io import fits
from scipy.stats import norm
from sklearn.neighbors import KernelDensity
import timeit

"""
################################################################################################################################
                    ################## SED fitting ########################
#################################################################################################################################

- Script to fit a modified blackbody profile to the 4 point SED generated using the residual radial profiles for the chosen source using MCMC. fitting method described in detail in Appendix 1 of the paper.  

###### Input files required ####################
1) .csv table files containing interpolated source residual profile and unc. data for each wavelength generated in the Radial Profiles step. eg:For CIT6 - cit6_res_interp.csv, cit6_res_interp_unc.csv
2) .dat table file containing the required radial points. x_interp.dat
3) .fits file of the beta profile of M31. Smith_etal_2012_M31_BetaMap.fits
4) .txt files containing filter response curves for SCUBA-2 and PACS - pacs_filter_response_curve_blue70.txt, pacs_filter_response_curve_red160.txt, scuba2_filter_response_curve_450.txt, scuba2_filter_response_curve_850.txt. 
	- SCUBA-2 filter response curves downloaded from the EAO JCMT SCUBA-2 filters page: "http://www.eaobservatory.org/jcmt/instrumentation/continuum/scuba-2/filters/". PACS filter response curves downloaded from the SVO Filter Profile Service (Rodrigo et al., 2012): "http://svo2.cab.inta-csic.es/svo/theory/fps/index.php?mode=browse&gname=Herschel&gname2=Pacs".
5) PACS PSFs in .fits form downloaded from Bocchio et al., 2016 via Vizier 

######### Output ####################

1) .csv file containing Temperature, dust mass column density and beta radial profile for the chosen source which can later be used to plot these profiles. 

"""

#Creating modified BB function (equation from Gordon et al., 2014)
def surface_brightness(frequencies, temperature, density, beta):
	
	
	#kappa = (9.6 / ((160e-6) ** (-1 * beta))) * ((c.value/frequencies) ** (-1 * beta)) #(9.6 / (c.vlaue/(160e-6) ** (-1 * beta)) - c.value taken out. //// * ((c.value/frequencies) ** (-1 * beta)) = E_lambda in Gordon et al., 2014. Gorden Original

	kappa = (kappa_eff / ((160e-6) ** (-1 * beta))) * ((c.value/frequencies) ** (-1 * beta)) #(9.6 / (c.vlaue/(160e-6) ** (-1 * beta)) - c.value taken out. //// * ((c.value/frequencies) ** (-1 * beta)) = E_lambda in Gordon et al., 2014. First kappa_eff=9.6 replaced with 26 which is more suitable kappa_eff for C-rich AGB stars like this one at 160micron. For O rich we need kappa_eff=8.8 at 160 micron. 
	
	S_nu = density * kappa * blackbody_nu(frequencies, temperature).to(u.Jy / (u.arcsec**2)).value #From Gordon et al., 2014. 
	
	return S_nu



#Building modifed BB (Must be global to be fed into multiple functions)
temperature = 40 * u.K #u.k = units in klevin
#wavelengths = np.array([70,160,450,850]) *u.um #once you attach units it becomes an astropy quantity - extension of an numpy array
wavelengths = 10**np.arange(1.,3,0.0005) *(u.um) #frequencies given in log values to make sure each increment is the same with increasing value. from 10micron(10^1) to 1000micron(10^3) in increments of 200(1/200).#IMPORTANT!!
frequencies = (c/wavelengths.to(u.m)).value #converting frequencies to Hz after it is converted from micron to m to meet m/s unit of c. #IMPORTANT!!
kappa_eff = 1 #Will be upadtaed by the global value at each iteration.


#Not really needed (only for testing purposes of code)
density = (-4) #(-4) = Logsclae value given. units = g/cm^2 - col.density (linear - 0.0001 g/cm^2). #Not really needed (only for testing purposes of code)
beta = 1.5 #1.5 #Not really needed (only for testing purposes of code)
bb = surface_brightness(frequencies, temperature, 10**density, beta) #Not really needed (only for testing purposes of code)
#bb = blackbody_nu(frequencies, temperature).to(u.Jy / (u.arcsec**2)) 
 



#Kernel Density Estimation of the beta profile from M31 in Smith et al., 2012 for the beta prior. (Must be global to be fed into multiple functions)
fitsFile = fits.open('Smith_etal_2012_M31_BetaMap.fits', memmap=True) #Importing fits file 
image = fitsFile[0].data 
image_numbers_only_array = image[~np.isnan(image)]

#Plotting data histogram
#nbins = 200
#histrogram = plt.hist(image_numbers_only_array, nbins, normed=1, color='blue') #normed=1: normalize the histogram to one (area of histogram=1) cause the KDE is auto normalised to one so in oder to view them in the same scale. 
#plt.show()

#Plotting progression to kernels
X = image_numbers_only_array[:,np.newaxis] #np.newaxis = Adds a new axis to the array. 1D array becomes 2D array, etc...
#print ('X =', X)

#Bandwidth 
n = len(X) # No. of data points 
d = 1 #Measuring only beta therefore 1D
bandwidth_calc = (n * (d+2) / 4.) ** (-1. / (d+4)) #Silverman's Rule

#Gaussina Kernel Density Estimation
#kde = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(X) 
#kde = KernelDensity(kernel='cosine', bandwidth=bandwidth_calc).fit(X) 
kde = KernelDensity(kernel='cosine', bandwidth=bandwidth_calc).fit(X) 
#kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth_calc*1.3).fit(X) #1.3 set by eye to find the best smoothing for the curve and minimize the histogram spikes as much as possible.
#log_dense = kde.score_samples(X_plot) #Insert into beta prior
#kde_plot = plt.plot(X_plot[:, 0], np.exp(log_dense), color='red')  




#Reading in fileter response curve files. 
def read_fileterresponse(frequencies):
	filenames = ("pacs_filter_response_curve_blue70.txt", "pacs_filter_response_curve_red160.txt", "scuba2_filter_response_curve_450.txt", "scuba2_filter_response_curve_850.txt") 

	#Convolving modified-bb with the filter response curves 
	filter_array = np.zeros([4,len(frequencies)]) #creating an empty 4 col(4 frequencies) array to hold each created synthetic photometry point once it is created 

	for i in range(len(filenames)):
		filename = filenames[i]

		if (i == 0) or (i == 1): #PACS filter response curves
			x_col, y_col = loadtxt(filename, usecols = (0,-1), unpack=True)	#x_col = col containing x axis values in filter response curve. y_col = col containing y axis values (dimensionless) in filter response curve. usecols = (1stcol, lastcol)
			x_col = (c/(x_col*1e-10)) #orginal units are in angstrom wavelength units. As our bb function is in Hz we convert this to Hz as well. As c is in m/s so x_col must be first converted from angstrom to m before converting to Hz. 

		else: #SCUBA2 filter response curves
			x_col, y_col = loadtxt(filename, usecols = (0,-1), unpack=True) #usecols = (1stcol, lastcol)
			x_col = (x_col*1e9) #orginal units are in GHz freq units. As our bb function is in Hz we convert this to Hz as well. 

		interpolated_function = interpolate.interp1d(x_col, y_col, kind='cubic', fill_value=0., bounds_error=False) #creating interpolating function using filter response data. 
		interpolated_filter = interpolated_function(frequencies) #applying interpolating function and interpolating to match the frequencies, i.e: x axis values of bb curve. In these two steps, we interpolate the filter response curves to match the x axis values for the bb curve so that the convolution later can work. 
		filter_array[i,:] = interpolated_filter#np.append(filter_array, interpolated_filter) #saving interpolated filter curves onto an array for later use. As the filter interpolated filter curves must be reused many times over it's easier to save it instead of running the loop over everytime. 
		
		#print interpolated_filter
       		#print np.sum(interpolated_filter > 0.01)
		#print interpolated_filter
		#print filter_array 
	
	return filter_array


#Reading in the filter array needed globally and setting it to a global parameter
filter_array = read_fileterresponse(frequencies) 



#Generating Synthetic Photometry
def synthetic_photometry(blackbody, filter_array):

	#for i in range(len(filenames)):
	syn_phot = np.sum(blackbody*filter_array, axis = 1) / np.sum(filter_array, axis = 1) #This is the SYNTHETIC PHOTOMETRY we generate for the model fitting. Here we convolve the bb spectrum with each filter response curve. "/np.sum..." = normalising the synthetic photometry
	#print syn_phot
	
	return syn_phot

	
def prior(theta): #Proiors - our data set. All in log values

	temperature, density, beta = theta
	
	#data for temperature normal distribution 	
	temperature_inner = 1300 #Temperature at the inner most radius of the star
	radius_inner = 0.03 #unit-arcsecond. Point at which the inner radius is located
	temperature_mean = (temperature_inner * (np.sqrt(radius_inner/x_point))) #x_point - in this case corresponds to the outermost radial point at that instance. The radial point which is being calculated at that instance. 
	temperature_sigma = 40 #uncertanity on temeprature #uncertainty in Kelvin. eg - 40 = +/-40 kelvin NOT 40%!
	
	#data for Beta from Global KDE 
	log_dense = kde.score_samples(beta)
		
	if 2.7 < temperature < 300 and -10 < density < 10 and 0 < beta:
		return norm.logpdf(temperature, temperature_mean, temperature_sigma) + log_dense
		
	return -np.inf 
	
	#OLD Priors 
	#if 2.7 < temperature < 800 and 0 < 10**density < np.inf and 0 < beta < 3:
	#if 2.7 < temperature < 300 and -10 < density < 10 and 0 < beta < 3:
	#	return 0.0
	#return -np.inf 	


#Building Model using Synthetic Photometry and bb curve equation.   
def likelihood_function(theta, y, yerr): #theta - the input values, x = frequencies, y = our fluxes from the residual profile. #all in log values

	temperature, density, beta = theta #input parameters for the bb function to make the MBB everytime over the chain
	model = surface_brightness(frequencies, temperature, 10**density, beta) #model in this instance = bb with the varying free parameters(frequencies, temperature, density, beta) given by mcmc
	model = synthetic_photometry(model, filter_array) #Model/Likleyhood function - Synthetic photometry is the model. We need to fit our data to this model.
 
	#We need two likelihood functions. One for detection data points (val>3sigma) and one for no dection data points(val<3sigma).
	#1. Gaussian function - Suitable for data with a well defined value with +/- unc and we have a detection of 3sigma or higher - value > 3*val.unc=3sigma.unc. Not good for values which are only upper limits.
	detection = np.where (y >= 3*yerr)
	likelihood = (-0.5) * np.sum ( ((y[detection]-model[detection])**2 / (yerr[detection] ** 2)) + (np.log(2*np.pi*(yerr[detection]**2))) ) 	
	#2. Cumilative Distribution function. For values where we only have an upper limit (In our data that's mostly 450 SCUBA2). To be used for data points where value < 3*value.unc=3sigma.unc. In this case 3yerr(3sgima value) which is the upper limit for that data point becomes the data point.
	non_detection = np.where (y < 3*yerr)	
	likelihood = likelihood + np.sum(norm.logcdf(3*yerr[non_detection], model[non_detection], yerr[non_detection]))

	return likelihood
	
	
def posterior(theta, y, yerr):

	prob = prior(theta)
	if not np.isfinite(prob):
		return -np.inf
	probablilty = prob + likelihood_function(theta, y, yerr)
	
	return probablilty




#Running MCMC and Saving the results
def run_emcee(flux_data, flux_unc_data, flux_dim, mcmcoutput_filename, ndim, nwalkers, steps, burn_in, radpoint):


	# Set up the Sampler
	pos=[(28.,1.,1.) + np.random.rand(ndim) for i in range(nwalkers)] #starting points/positions for temperature(28K + some random perturbation), density(1g/cm2 + some random perturbation), beta(1 + some random perturbation) walkers in the moddle grid. This number doesn't really matter as walkers will then move around and converge on the best correct number. ndim&nawalkers given in main body.
	#sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior, args=(synthetic_photometry(bb, filter_array), 0.1*synthetic_photometry(bb, filter_array))) #Here we use synthetic photometry data as our sample set of data and unc= 10%*synthetic photometry values.
	sampler = emcee.EnsembleSampler(nwalkers, ndim, posterior, args=(flux_data, flux_unc_data))#Here we use the flux data from our radial profiles as our sample set of data and unc.


	# Run the production chain.
	print("Running MCMC...") 
	sampler.run_mcmc(pos, steps) #steps given in main body.
	print("Done.")
	 


	#Plotting Walker plots
	fig, axes = plt.subplots(3, 1, sharex=True, figsize=(8, 9))
	
	axes[0].plot(sampler.chain[:, :, 0].T, color="k", alpha=0.4)
	axes[0].set_ylabel("$temp$")

	axes[1].plot(sampler.chain[:, :, 1].T, color="k", alpha=0.4)
	axes[1].set_ylabel("$10**density$")

	axes[2].plot((sampler.chain[:, :, 2]).T, color="k", alpha=0.4)
	axes[2].set_ylabel("$beta$")
	axes[2].set_xlabel("step number")
 
	plt.savefig(radpoint + '_walkerplot.png')	
	#plt.show()


	#Saving the data for Plotting corner plots and analysis
	samples = sampler.chain[:, burn_in:, :].reshape((-1, ndim)) #sampler.chain[no.of.chains, burn-in-step, no.of.variables]. #burnin - given in main body - point/step from which we are to consider the data. The point where the walkers all converge as seen in the walker plots.


	#Save the data after the burnin has taken place. 
	f = open(mcmcoutput_filename, "w") #creating and wrtting(w) to a file. mcmcoutput_filename - given in mainbody. 
	f.close()

	f = open(mcmcoutput_filename, "a") #append(a) the created data file 
	for i in range(samples.shape[0]):
		outstring = "\t" + str(samples[i,0]) + "\t" + str(samples[i,1]) + "\t" + str(samples[i,2]) + "\t" + "\n" #"\t" = tab space. "\n" = start new line. str() = converting/casting a value into a string-Textual data in Python is handled with str objects. We're just printing the data as text format(allthough they are numbers) into the file.
		f.write(outstring)
	f.close()

	#return - In this case since there is already an output generated we do not need a "return".


#Analysing emcee results and plotting corner plots
def analyse_emcee(mcmcoutput_filename, ndim, radpoint): 

	
	data = np.loadtxt(mcmcoutput_filename) #Loading the file contaning the temp, density, beta array data from mcmc
	samples = data.reshape([len(data),ndim]) #Reshaping the saved array to a '3xlength' array (in this case it's the same shape as the orginal array) that can be read and used to plot corner plots. #ndim - Given in main body. len(data) - length of the data array. 


	#Plotting Corner plots
	fig = corner.corner(samples, labels=["$temperature$", "$log(density)$", "$beta$"])
	#fig = corner.corner(samples, labels=["$temperature$", "$log(density)$", "$beta$"], truths=[temperature.value, density, beta]) #truths : iterable (ndim,) A list of reference values to indicate on the plots.  Individual values can be omitted by using ``None``. truth_color : str A ``matplotlib`` style color for the ``truths`` makers.
	#fig = corner.corner(samples, levels=[0,0], labels=["$temperature$", "$log(density)$", "$beta$"], truths=[temperature.value, density, beta]) #levels=[0,0] - No contour lines
		
	fig.savefig(radpoint + '_cornerplot.png')
	#plt.show()


	#Final values and Uncertainties for temperature, density and beta
	temperature_mcmc, density_mcmc, beta_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0))) #*np.percentile(samples, [16-v0, 50-v1, 84-v2] = symmetrical/nornal distribution arround the median assuming the distribution is gaussian. 50th percentile - median=our result. 84th percentile - +1sigma value(val at 84v2 - val at 50v1 = our +unc.). 16th percentile - -1sigma value(val at 50v1 - val at 16v0 = our -unc.). Look at 3sigma plot in https://en.wikipedia.org/wiki/Percentile. 
	
	#print temperature_mcmc
	#print density_mcmc
	#print beta_mcmc
	
	
	return  temperature_mcmc, density_mcmc, beta_mcmc
	

#Saving the output list into an array (so it can be saved in a way that can be reread for plotting and then the array into a file.	
def saving_mcmc_output(temperature_save, density_save, beta_save, temperature_filename, density_filename, beta_filename, x_interp):

	temperature_array = np.asarray(temperature_save) #Casting a list into an array. list-temperature_save. command-'np.asarray'
 	temperature_reshaped = temperature_array.reshape((len(temperature_array),3)) #Reshaping the array to a legth_of_array x 3 (3D=val,+unc,-unc.)
 		
 	f = open(temperature_filename, "w") 
	f.write("x_interp" + "," +  "Temperature" + "," + "Plus_Unc" + "," + "Minus_Unc" + "\n")
	f.close()

	f = open(temperature_filename, "a") 
	for i in range(temperature_reshaped.shape[0]):
		outstring = str(np.float64(x_interp[i])) + "," + str(temperature_reshaped[i,0]) + "," + str(temperature_reshaped[i,1]) + "," + str(temperature_reshaped[i,2]) + "\n" 
		f.write(outstring)	
	f.close()

 	density_array = np.asarray(density_save) #Casting a list into an array. list-densitysave. command-'np.asarray'
 	density_reshaped = density_array.reshape((len(density_array),3)) #Reshaping the array to a legth_of_array x 3 (3D=val,+unc,-unc.)
 		
 	f = open(density_filename, "w") 
	f.write("x_interp" + "," + "Density" + "," + "Plus_Unc" + "," + "Minus_Unc" + "\n")
	f.close()

	f = open(density_filename, "a") 
	for i in range(density_reshaped.shape[0]):
		outstring = str(np.float64(x_interp[i])) + "," + str(density_reshaped[i,0]) + "," + str(density_reshaped[i,1]) + "," + str(density_reshaped[i,2]) + "\n" 
		f.write(outstring)	
	f.close()

	beta_array = np.asarray(beta_save) #Casting a list into an array. list-beta_save. command-'np.asarray'
 	beta_reshaped = beta_array.reshape((len(beta_array),3)) #Reshaping the array to a legth_of_array x 3 (3D=val,+unc,-unc.)
 		
 	f = open(beta_filename, "w") 
	f.write("x_interp" + "," + "Beta" + "," + "Plus_Unc" + "," + "Minus_Unc" + "\n")
	f.close()

	f = open(beta_filename, "a") 
	for i in range(beta_reshaped.shape[0]):
		outstring = str(np.float64(x_interp[i])) + "," + str(beta_reshaped[i,0]) + "," + str(beta_reshaped[i,1]) + "," + str(beta_reshaped[i,2]) + "\n" 
		f.write(outstring)	
	f.close()
	
	#return temperature_filename, density_filename, beta_filename
	

#Loading x axis from the x_interp.dat file
def radial_point(star):

	x_interp_data = np.loadtxt('x_interp.dat')
	#x_interp_data = np.loadtxt(star + '_x_interp_test.dat') # x_interp_test file only for testing DOES NOT HAVE THE FULL X AXIS ARRAY!!!!!!
	x_interp = x_interp_data.reshape([len(x_interp_data),1]) #(len,col) 
	

	#print x_interp
	#sys.exit("Forced Exit by user")
	
	return x_interp


	

#Plotting Temperature, Density, Beta Radial Profiles. 
def temp_dens_beta_radial_profiles(star, temperature_filename, density_filename, beta_filename, nskip=1): #nskip=1 - skip the first radial point (0th) and plot the rest

	
	#x_interp = radial_point(star) #Loading x_interp values(x radial points) from the function which calls on the file with the saved x_axis values. 
	#x_interp = x_interp[1:] #x_interp is now a 1D array and we take all values from the first(2nd) point, NOT the zeroth(1st) point.
	

	temperature_data = Table.read(temperature_filename, format='ascii.csv') #Loading the file contaning the temp, density, beta array data from mcmc
	temperature_x_interp = temperature_data['x_interp']
	temperature_y = temperature_data['Temperature'] #Temperature values
	temperature_plusunc =  temperature_data['Plus_Unc'] # + uncertainties
	temperature_minusunc = temperature_data['Minus_Unc'] # - uncertainties


	density_data = Table.read(density_filename, format='ascii.csv') 
	density_x_interp = density_data['x_interp']
	density_y = density_data['Density'] #density values
	density_plusunc =  density_data['Plus_Unc'] # + uncertainties
	density_minusunc = density_data['Minus_Unc'] # - uncertainties

	beta_data = Table.read(beta_filename, format='ascii.csv')
	beta_x_interp = beta_data['x_interp']
	beta_y = beta_data['Beta'] #beta values
	beta_plusunc =  beta_data['Plus_Unc'] # + uncertainties
	beta_minusunc = beta_data['Minus_Unc'] # - uncertainties


	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.errorbar(temperature_x_interp[nskip:], temperature_y[nskip:], yerr=[temperature_minusunc[nskip:], temperature_plusunc[nskip:]], fmt='--o', color='orange')
	ax.set_xlim([0,150])
	#ax.set_ylim([0,250])
	ax.set_xlabel("Radius ($^{\prime\prime}$)")
	ax.set_ylabel('Temperature (K)')
	plt.title(star + "\t" + 'Temperature Profile')
	plt.savefig(star+'_TemperatureProfile.png')

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.errorbar(density_x_interp[nskip:], density_y[nskip:], yerr=[density_minusunc[nskip:], density_plusunc[nskip:]], fmt='--o', color='purple')
	ax.set_xlim([0,150])
	#ax.set_ylim([-10, -4])
	ax.set_xlabel("Radius ($^{\prime\prime}$)")
	ax.set_ylabel('log(col.density) (g/cm**2)')
	plt.title(star + "\t" + 'Density Profile')
	plt.savefig(star+'_DensityProfile.png')

	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.errorbar(beta_x_interp[nskip:], beta_y[nskip:], yerr=[beta_minusunc[nskip:], beta_plusunc[nskip:]], fmt='--o', color='green')
	ax.set_xlim([0,150])
	#ax.set_ylim([0,3])
	ax.set_xlabel("Radius ($^{\prime\prime}$)")
	ax.set_ylabel('beta')
	plt.title(star + "\t" + 'Beta Profile')
	plt.savefig(star+'_BetaProfile.png')

	#plt.show()
	
	return x_interp



x_point = 0 #x_point = radial point at which the calculation is carried out. Set to zero for now and will be overwritten by the given radial point in the mcmc loops.
 



################### Main Body - Script ##############################################


if __name__=="__main__":

	star = 'uhya'
	global kappa_eff
	kappa_eff = 26 #Kappa_eff = 26 for C rich sources and Kappa_eff = 8.8 for O rich and S-type and RSGs
	flux_file = star + '_res_interp.dat'
	flux_unc_file = star + '_res_interp_unc.dat'
	flux_dim = 4
	ndim = 3 # ndim=no. of dimensions
	nwalkers = 100 #nwalkers=no. of walkers (Fast Run: nwalkers=30, steps=1000, burn_in=200 // Long Run: nwalkers=100, steps=10000, burn_in=1000)
	steps = 10000
	burn_in = 1000 #burnin - point/step from which we are to consider the data. The point where the walkers all converge as seen in the walker plots. 
	temperature_filename = star + '_temperature_radial_output.csv'
	density_filename = star + '_density_radial_output.csv'
	beta_filename = star + '_beta_radial_output.csv'
	

	
	
	#Opening the files containing Flux data from Radial profiles for the SED fitting
	flux_data_load = np.loadtxt(flux_file)
	flux_data = flux_data_load.reshape([len(flux_data_load), flux_dim])

	flux_unc_data_load = np.loadtxt(flux_unc_file)
	flux_unc_data = flux_unc_data_load.reshape([len(flux_unc_data_load), flux_dim])
	
	x_interp = radial_point(star)
	
	for i in range(0, len(flux_data)): #0 = start from 1st(0) row.
	#for i in range(0, 5): #Testing purposes only!! For all radial points use the other one
	
	
		radpoint = str(i)
		mcmcoutput_filename = radpoint + '_emcee.dat'		
		x_point = x_interp[i] #Outer radial point
		
		
		
		
		start_time = timeit.default_timer() #Start timer 
		print ('start time =', start_time)
		run_emcee(flux_data[i,:], flux_unc_data[i,:], flux_dim, mcmcoutput_filename, ndim, nwalkers, steps, burn_in, radpoint)
		temperature_output, density_output, beta_output = analyse_emcee(mcmcoutput_filename, ndim, radpoint) #We put the outputs(returns) from the function 'analyse_emcee' into the variables listed. The retruns will be put out in the order they are given in the 'return' command in the function, so the variables have to be ordered accordingly. 
		elapsed = timeit.default_timer() - start_time #End timer
		print ('end time =', elapsed) 
		

		#Saving the outputs into lists so they can be used for radial profile plotting. This bit is in the loop because the list has to be appended each time around the loop. 
		try: 
 			temperature_save.append([temperature_output])
 		except:
 			temperature_save = ([])
 			temperature_save.append([temperature_output])

		try: 
 			density_save.append([density_output])
 		except:
 			density_save = ([])
 			density_save.append([density_output])

		try: 
 			beta_save.append([beta_output])
 		except:
 			beta_save = ([])
 			beta_save.append([beta_output])
 			
 		
 	
 	#Saving the output list into an array (so it can be saved in a way that can be reread for plotting and then the array into a file.
 	saving_mcmc_output(temperature_save, density_save, beta_save, temperature_filename, density_filename, beta_filename, x_interp)
 	
 	#Plotting Temperature, Density, Beta Radial Profiles.
 	temp_dens_beta_radial_profiles(star, temperature_filename, density_filename, beta_filename, nskip=1)




	 


 


















