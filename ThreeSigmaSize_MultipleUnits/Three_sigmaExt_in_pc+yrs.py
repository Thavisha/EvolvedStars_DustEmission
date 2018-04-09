import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import pylatex
import matplotlib.ticker as mtick

"""
###################################################################################################################################
    ########### Three Sigma Radii in Multiple Units  ################################
###################################################################################################################################

- Convert the Three sigma extents for all four wavelengths measured from the Residual profiles into physical size (pc) and age of CSE extenet (yr)


###### Input files required ####################

1) Source_Info.csv containing the details of U Ant. 

##### Output ############

.csv table containing the input source details and the output in multiple units for all three wavelengths. 

"""

Source_Data = Table.read('Source_Info.csv', format='ascii.csv')

f = open('threeSigma_Radius_Multiple_Units.csv', 'w')
f.write("source" + "," + "distance" + "," + "Terminal_Velocity_cm/s" + "," + "PACS70_threeSigma_arcsec" + "," + "PACS160_threeSigma_arcsec" + "," + "SCUBA450_threeSigma_arcsec" + "," + "SCUBA850_threeSigma_arcsec" + "," + "PACS70_threeSigma_pc" + "," + "PACS160_threeSigma_pc" + "," + "SCUBA450_threeSigma_pc" + "," + "SCUBA850_threeSigma_pc" + "," + "PACS70_threeSigma_yrs" + "," + "PACS160_threeSigma_yrs" + "," + "SCUBA450_threeSigma_yrs" + "," + "SCUBA850_threeSigma_yrs" + "," + "\n")
f.close()

for Source in Source_Data:
	
	star = Source['Source'] 			
	distance = Source['Distance(pc)'] #units pc
	velocity = Source['DeBeck_etal_2010_Terminal_Velocity_fromCOdata(cm/s)']  #Terminal/Outflow velocity || units=cm/s || From DeBeck et al., 2010
	seventy_3sigma_arcsec = Source['PACS70_threeSigma_arcsec'] #PACS 70 three sigma in arcsec units
	onesixty_3sigma_arcsec = Source['PACS160_threeSigma_arcsec'] #PACS 160 three sigma in arcsec units
	fourfifty_3sigma_arcsec = Source['SCUBA450_threeSigma_arcsec'] #SCUBA-2 450 three sigma in arcsec units
	eightfifty_3sigma_arcsec = Source['SCUBA850_threeSigma_arcsec'] #SCUBA-2 850 three sigma in arcsec units
	
	#Three sigma radius from arcsec (projected) to pc (physical)
	seventy_3sigma_pc = (seventy_3sigma_arcsec * distance) * 4.84814e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc
	onesixty_3sigma_pc = (onesixty_3sigma_arcsec * distance) * 4.84814e-6
	fourfifty_3sigma_pc = (fourfifty_3sigma_arcsec * distance) * 4.84814e-6
	eightfifty_3sigma_pc = (eightfifty_3sigma_arcsec * distance) * 4.84814e-6

	
	#Three sigma radius from pc to cm for age calc. 
	seventy_3sigma_cm = seventy_3sigma_pc * 3.08567961169e+18 
	onesixty_3sigma_cm = onesixty_3sigma_pc * 3.08567961169e+18 
	fourfifty_3sigma_cm = fourfifty_3sigma_pc * 3.08567961169e+18 
	eightfifty_3sigma_cm = eightfifty_3sigma_pc * 3.08567961169e+18 

	#Three sigma radius in time (yrs) units for age. 
	seventy_3sigma_yrs = (seventy_3sigma_cm / velocity) / 3.154e+7 #/3.154e+7 to convert from s to years
	onesixty_3sigma_yrs = (onesixty_3sigma_cm / velocity) / 3.154e+7 
	fourfifty_3sigma_yrs = (fourfifty_3sigma_cm / velocity) / 3.154e+7 
	eightfifty_3sigma_yrs = (eightfifty_3sigma_cm / velocity) / 3.154e+7 


	f = open('threeSigma_Radius_Multiple_Units.csv', 'a')
	outstring = str(star)  + ","  + str(distance)  + ","  + str(velocity)   + ","  + str(seventy_3sigma_arcsec) + "," + str(onesixty_3sigma_arcsec) + "," + str(fourfifty_3sigma_arcsec) + "," + str(eightfifty_3sigma_arcsec) + "," + str(seventy_3sigma_pc) + "," + str(onesixty_3sigma_pc) + "," + str(fourfifty_3sigma_pc) + "," + str(eightfifty_3sigma_pc) + "," + str(seventy_3sigma_yrs) + "," + str(onesixty_3sigma_yrs) + "," + str(fourfifty_3sigma_yrs) + "," + str(eightfifty_3sigma_yrs) + "\n"
	f.write(outstring)	
	f.close()


	#seventy_arcsec_median = np.median(seventy_3sigma_arcsec)
	#seventy_arcsec_25thpercentile = np.percentile(seventy_3sigma_arcsec, 25)
	#seventy_arcsec_75thpercentile = np.percentile(seventy_3sigma_arcsec, 75)
	#seventy_iqrange = seventy_arcsec_75thpercentile - seventy_arcsec_25thpercentile

	#onesixty_arcsec_median = np.median(seventy_3sigma_arcsec)
	#onesixty_arcsec_25thpercentile = np.percentile(seventy_3sigma_arcsec, 25)
	#onesixty_arcsec_75thpercentile = np.percentile(seventy_3sigma_arcsec, 75)
	#onesixty_iqrange = seventy_arcsec_75thpercentile - seventy_arcsec_25thpercentile









































