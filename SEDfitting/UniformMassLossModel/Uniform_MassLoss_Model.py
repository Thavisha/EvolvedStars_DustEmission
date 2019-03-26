import matplotlib.pyplot as plt
import numpy as np
from numpy import trapz
from scipy.integrate import simps
from astropy.table import Table

"""
###################################################################################################################################
    ########### Uniform Mass Loss model ################################
###################################################################################################################################

- Generating a radially dependent dust mass column density profile assuming uniform mass loss. 

- We assume dust mass column density \prop to 1/r^2 relationship and extending it in the z and impact parameter directions. 


###### Input files required ####################

 - x_interp.dat file containing x axis information for upto which to generate the UMLM

###### Output ##################################  

1) .csv file containing the dust mass column density radial profile derived assuming uniform mass loss. This profile can now be over-plotted with the dust mass column density derived using SED fitting to show the deviation of the SED derived profile from uniform mass loss. 

"""

def Uniform_MassLoss_Model(impact_para, impact_para_0, z_0, r_max): 
#Project radius vs. Density is what we get out. Using density is proportional to 1/r^2 law and extending it in the z and impact parameter directions,

	impact_para = impact_para

	z_max = np.sqrt((r_max ** 2) - (impact_para ** 2))

	z = np.arange(0,z_max)

	density = (
			(np.sqrt((impact_para_0 ** 2) + (z_0 ** 2)))

				/(np.sqrt((impact_para ** 2) + (z ** 2)))

			) ** 2

	return density, z

f = open('Uniform_MassLoss_Model.csv', "w") #All other sources except IRC+10216
#f = open('Uniform_MassLoss_Model_longXaxisforIRC10216.csv', "w") #Long x axis for IRC+10216
f.write("x_interp" + "," + "Uniform_MassLoss_Model_Values" + "\n") 
f.close()



if __name__=="__main__":

	surface_density = [] #Defining an empty list to save the surface density output

	x_interp_data = np.loadtxt('x_interp_Uniform_MassLoss_Model.dat') # All sources except IRC+10216
	#x_interp_data = np.loadtxt('x_interp_irc10216_longxaxis_Uniform_MassLoss_Model.dat') # Long x axis for IRC+10216
	x_interp = x_interp_data.reshape([len(x_interp_data),1]) 

	for value in x_interp[:-1]: 
		impact_para = value

		impact_para_0 = 20
		z_0 = impact_para_0

		r_max = x_interp[-1]

		density, z = Uniform_MassLoss_Model(impact_para, impact_para_0, z_0, r_max)				
		
		surface_density.append(simps(density, z)) #Calculating surface density and saving it onto the blank surface_density list
	
	f = open('Uniform_MassLoss_Model.csv', "a") #All other sources except IRC+10216
	#f = open('Uniform_MassLoss_Model_longXaxisforIRC10216.csv', "a") #Long x axis for IRC+10216
	for i in range(len(surface_density)):
		outstring = str(x_interp[i,0]) + "," + str(surface_density[i]) + "\n" 
		f.write(outstring)
	f.close()
	
	#print x_interp[:,0].shape
	#print x_interp[:,0], surface_density,x_interp[1,0]

	

	#Plotting density and mass loss model together for stellar source

	star_density_table = Table.read('cit6_density_radial_output.csv', format='ascii.csv')
	star_density = star_density_table['Density']
	star_density_plus_unc =  star_density_table['Plus_Unc']
	star_density_minus_unc =  star_density_table['Minus_Unc']
	
	#print np.array(surface_density)*surface_density[2]/10**star_density[2], star_density.data[2]

	surface_density_scaled = ((10.**(star_density[2]))/surface_density[3])*np.array(surface_density) #Scaling uniform mass loss model to match a point on the stellar density profile. In this case we are matching(scalling to) the third(2) radial points. #Stellar density converted from log to linear for scaling

	fig = plt.figure(figsize=(8, 9)) #(width,height)
	gs = gridspec.GridSpec(1, 1) #Note capital G and S in .GridSpec
		
	ax1 = fig.add_subplot(111)#gs[0]) 
	ax1.errorbar(x_interp[1:-1], star_density, yerr=[star_density_minus_unc, star_density_plus_unc], fmt='--^', color='indigo')
	ax1.plot(x_interp[:-1], np.log10(surface_density_scaled), '--', color='red') #Plotting surface density model in log form
	ax1.set_ylim([-10, -1])
	ax1.set_ylabel('log($\\Sigma$ (g cm$^{-2}$))')
	ax1.set_xlabel('Radius (")')

	#plt.savefig("CIT6_Density+UniformMassLossModel_test1.png")
	#plt.show()

	




