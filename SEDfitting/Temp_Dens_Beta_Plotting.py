import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import pylatex
import matplotlib.ticker as mtick

"""
################################################################################################################################
                    ################## Plotting resulting Temperature, Density and Beta Profiles ########################
#################################################################################################################################

- script which plots the radially dependent temperature, density and beta profiles resulting from the SED mcmc fitting.
- The profile x axes are presented in both projected radius (arcsec) and age of CSE (years). 
- The age is calculated using the physical radius (cm - derived using projected radius and distance to source) and the terminal velocity (from De Beck et al., 2010) - (time = dist / vel) 

###### Input files required ####################
1) .csv file containing source names, distances, terminal velocities, Anchor points for uniform mass loss model - Source_Information.csv
2) .csv table file containing the required radial points. x_interp.csv
3) .csv file containing Temperature, dust mass column density and beta radial profiles for all sources from the above. eg: for CIT6 - cit6_temperature_radial_output.csv, cit6_density_radial_output.csv, cit6_beta_radial_output.csv.
4) .csv file containing the dust mass column density radial profile derived assuming uniform mass loss derived from above: uniform_MassLoss_Model.csv - for over-plotting on the mcmc derived  dust mass column density to show the deviation of it from uniform mass loss. 
	- This profile is scaled to a chosen point on the  mcmc derived  dust mass column density in order to plot it.   

######### Output ####################
1) .png of image of the temperature, density and beta profile for all the sources in the sample. 


"""

font = {'family' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

#Change depending on Source!!!!! 	
star = 'whya' #	
distance = 104.275  #units pc #			
velocity = 850000 #Terminal/Outflow velocity || units=cm/s || From DeBeck et al., 2010
radial_point_to_anchor_UMLM = 3 #Radial point at which to anchor Uniform Mass Loss Model (UMLM) on the density profile. Counting starts from 0 in python for the radial points. So 1st radial point = 0, 2nd radial point = 1, etc.... 
radial_point_to_anchor_UMLM = radial_point_to_anchor_UMLM + 1 #Since we're skipping plotting the 0th radial point so the UMLM needs to be anchord to the 2nd plotted radial point. 
nskip = 1 #Skip first radial point (0th) in all profiles and plot from the second(1) radial point.  


#Loading and adjusting x_interp
#x_interp in arcsec
x_interp_data = Table.read('x_interp.csv', format='ascii.csv')
x_interp = x_interp_data['x_interp'] #Units in arcsec

#x_interp in cm 
x_interp_cm = (distance * x_interp) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm

#x_interp in time (s) units (time = dist / vel) 
x_interp_time = x_interp_cm / velocity #units in s
x_interp_time = x_interp_time / 3.154e+7 #units in years (1yr = 3.154e+7)


#Loading Temperature, Density and Beta Profiles
temperature_filename = star + '_temperature_radial_output.csv'
temperature_data = Table.read(temperature_filename, format='ascii.csv')
Temperature = temperature_data['Temperature']
Temperature_PlusUnc = temperature_data['Plus_Unc']
Temperature_MinusUnc = temperature_data['Minus_Unc']

density_filename = star + '_density_radial_output.csv'
density_data = Table.read(density_filename, format='ascii.csv')
Density = density_data['Density']
Density_PlusUnc = density_data['Plus_Unc']
Density_MinusUnc = density_data['Minus_Unc']

beta_filename = star + '_beta_radial_output.csv'
beta_data = Table.read(beta_filename, format='ascii.csv')
Beta = beta_data['Beta']
Beta_PlusUnc = beta_data['Plus_Unc']
Beta_MinusUnc = beta_data['Minus_Unc']


#Loading Uniform Mass Loss Model (UMLM) and scaling it (#0th radial point(=inf) deleted from Model to make it easy for plotting)
Uniform_MassLoss_Model_file = Table.read('Uniform_MassLoss_Model.csv', format='ascii.csv')
Uniform_MassLoss_Model = Uniform_MassLoss_Model_file['Uniform_MassLoss_Model_Values']

#Scaling Uniform Mass Loss Model to anchor (match) at a given radial point 
Uniform_MassLoss_Model_scaled = ((10.**(Density[radial_point_to_anchor_UMLM]))/Uniform_MassLoss_Model[radial_point_to_anchor_UMLM])*np.array(Uniform_MassLoss_Model) #Scaling uniform mass loss model to match a point [] on the stellar density profile. #Stellar density converted from log to linear for scaling


#Plotting
fig = plt.figure(figsize=(8, 9)) #(width,height)
gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 1]) #Note capital G and S in .GridSpec
gs.update(hspace=0)

ax1 = fig.add_subplot(gs[0]) #since we only have 1 col, we have to only give the row into which the plot must go
ax1.errorbar(x_interp[nskip:], Temperature[nskip:], yerr=[Temperature_MinusUnc[nskip:], Temperature_PlusUnc[nskip:]], fmt='--o', color='crimson', capsize=2)
nbins = len(ax1.get_yticklabels())
ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
#ax1.get_yaxis().set_label_coords(-0.01,0.5)
plt.setp(ax1.get_xticklabels(), visible=False)
temperature_lower_limit = (min(Temperature) - Temperature_MinusUnc[np.argmin(Temperature)]) - 20
temperature_upper_limit = (max(Temperature) + Temperature_PlusUnc[np.argmax(Temperature)]) + 20
#ax1.set_ylim([temperature_lower_limit, temperature_upper_limit])
#ax1.set_ylim([0,250])
#ax1.set_xlabel("Radius ($^{\prime\prime}$)")
ax1.set_ylabel('T (K)')


ax2 = fig.add_subplot(gs[1], sharex=ax1) 
ax2.errorbar(x_interp[nskip:], Density[nskip:], yerr=[Density_MinusUnc[nskip:], Density_PlusUnc[nskip:]], fmt='--^', color='indigo', capsize=2)
ax2.plot(x_interp[nskip:], np.log10(Uniform_MassLoss_Model_scaled[nskip:]), '--', color='orange', lw=3) #Plotting the scaled UMLM in log form
plt.setp(ax2.get_xticklabels(), visible=False)
nbins = len(ax2.get_yticklabels())
ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1f'))
ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis
#ax2.set_ylim(ymin=0) #starting y axis from zero
density_lower_limit = (min(Density) - Density_MinusUnc[np.argmin(Density)]) - 0.4
density_upper_limit = (max(Density) + Density_PlusUnc[np.argmax(Density)]) + 0.6
#ax2.set_ylim([density_lower_limit, density_upper_limit])
#ax2.set_ylim([-10, -4.5])
ax2.set_ylabel('log($\\Sigma$ (g cm$^{-2}$))')

	
ax3 = fig.add_subplot(gs[2], sharex=ax1) 
ax3.errorbar(x_interp[nskip:], Beta[nskip:], yerr=[Beta_MinusUnc[nskip:], Beta_PlusUnc[nskip:]], fmt='--s', color='darkgreen', capsize=2)
#plt.setp(ax3.get_xticklabels(), visible=False)
nbins = len(ax3.get_yticklabels())
ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis
nbins = len(ax3.get_xticklabels())
ax3.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis
ax3.set_xlim([0,145])
beta_lower_limit = (min(Beta) - Beta_MinusUnc[np.argmin(Beta)])
beta_upper_limit = (max(Beta) + Beta_PlusUnc[np.argmax(Beta)]) + 0.4
ax3.set_ylim([beta_lower_limit, beta_upper_limit])
#ax3.set_ylim([0,4])
ax3.set_ylim(ymin=0) #starting y axis from zero
ax3.set_xlabel("Radius ($^{\prime\prime}$)")
ax3.set_ylabel('$\\beta$')

#Additional x axis for time units
ax4 = ax1.twiny()
ax4.errorbar(x_interp_time[nskip:], Temperature[nskip:], color='none', mec='none')
temperature_lower_limit_2 = (min(Temperature) - Temperature_MinusUnc[np.argmin(Temperature)]) - 20
temperature_upper_limit_2 = (max(Temperature) + Temperature_PlusUnc[np.argmax(Temperature)]) + 20
#ax4.set_ylim([temperature_lower_limit_2, temperature_upper_limit_2])
#ax4.set_ylim([0,250])
#ax4.set_xlim([0,3.3])
nbins = len(ax4.get_xticklabels())
ax4.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on y axis to prevent merging with upper axis
#ax4.set_xlabel("Time ($\\times 10^{11} $s)")
#ax4.set_xlabel("Time ($10^{3}$yr)")
ax4.set_xlabel("Time (yr)")



save_file = star + '_Temp_Dens_Beta_timeXaxis+UMLM.png'
plt.savefig(save_file)
plt.show()




