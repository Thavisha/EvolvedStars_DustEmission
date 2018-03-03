import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
import scipy.optimize as optimization

"""
################################################################################################################################
################## Plotting Gas MLR vs. Dust MLR - with best fit line and Dust-to-Gas ratio = 1/200 line ########################
#################################################################################################################################

- Script to plot the CO Gas MLR vs. Dust MLR (dervied from integrating the density profile) with Chi sqaure minimized line of best fit and the canonical dust-to-gas ratio of 1/200 line. 

- Plotting results derived from the script described above. 

###### Input files required ####################
1) GasMLR_vs_DustMLR_withAccepted+BestFitLines.csv - table containing the information required to produce the plot. 


########## Output ##############
- plot of  CO Gas MLR vs. Dust MLR with Chi sqaure minimized line of best fit and the canonical dust-to-gas ratios for C and O rich.
"""


font = {'family' : 'normal',
        'size'   : 15}

plt.rc('font', **font)

#Defining function for fitting a chi sqaured minimization for the dust MLR vs. Gas MLR plot to show the best fit line. 
def line_function(x, b, m): #log(y)=mx+b therefore y=bx**m || x values has to be the first argument for curve fit!!
	function = b*(x**m)
	return function


#Reading in Source Data table 
source_data = Table.read('Plot_Data.csv', format='ascii.csv')

Source = source_data['source']

Stellar_Type = source_data['stellar_type']
o_rich = Stellar_Type == 'o' #Source[Stellar_Type == 'o']
c_rich = Stellar_Type == 'c'
rsg_type = Stellar_Type == 'rsg'
s_type = Stellar_Type == 's'

PACS160_three_sigma_arcsec = source_data['160_three_sigma(arcsec)']
PACS160_three_sigma_age_years = source_data['3sigma_radius_160_Time(years)']
Gas_MassLoss_Rate = source_data['DeBeck_CO_Mass_Loss_Rate(M_sun/yr)']
Dust_MassLoss_Rate = source_data['Dust_Mass_Loss_Rate(M_sun/yr)']
Dust_to_Gas_Ratio = source_data['Dust-to-gasRatio']


#Gas MLR vs. Dust MLR - with best fit line and Dust-to-Gas ratio = 1/200 line
mean_dust_to_gas_ratio = np.average(Dust_to_Gas_Ratio)
best_fit, covariance = optimization.curve_fit(line_function, Gas_MassLoss_Rate, Dust_MassLoss_Rate, [mean_dust_to_gas_ratio,1]) 
best_fit_y = (best_fit[0]*Gas_MassLoss_Rate)**best_fit[1] #b(x**m). b is the 0th element in best fit and m is the 1st. 
#print best_fit #[b=0.46310733,  m=1.31520542]
fig = plt.figure(figsize=(8, 6)) #(width,height)
plt.plot(Gas_MassLoss_Rate[o_rich], Dust_MassLoss_Rate[o_rich], 'o', markersize=10, color='mediumblue')
plt.plot(Gas_MassLoss_Rate[c_rich], Dust_MassLoss_Rate[c_rich], '*', markersize=10, color='crimson')
plt.plot(Gas_MassLoss_Rate[rsg_type], Dust_MassLoss_Rate[rsg_type], '^', markersize=10, color='maroon')
plt.plot(Gas_MassLoss_Rate[s_type], Dust_MassLoss_Rate[s_type], 's', markersize=10, color='forestgreen')
plt.plot(np.sort(Gas_MassLoss_Rate), 0.003*np.sort(Gas_MassLoss_Rate), '--', markersize=6, color='orange') #1/400 (=0.003) line. y=0.003x - C-rich
plt.plot(np.sort(Gas_MassLoss_Rate), 0.007*np.sort(Gas_MassLoss_Rate), '--', markersize=6, color='navy') #1/160 (=0.007) line. y=0.007x - O-rich
#plt.plot(np.sort(Gas_MassLoss_Rate), np.sort(best_fit_y), '--', markersize=6, color='red') #Line of best fit
plt.yscale('log')
plt.xscale('log')
plt.ylabel("Dust Mass Loss Rate ($M_{\\odot}/yr$)")
plt.xlabel("CO Mass Loss Rate ($M_{\\odot}/yr$)")
plt.savefig('GasMassLossRate_Vs_DustMassLossRate_C-ORatioLines.png')
#plt.savefig('GasMassLossRate_Vs_DustMassLossRate_withBestFitLin.png')
plt.show()




























