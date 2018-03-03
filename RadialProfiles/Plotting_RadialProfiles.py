import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import pylatex
import matplotlib.ticker as mtick
import warnings

"""
################################################################################################################################
                    ################## Plotting Radial Profiles ########################
#################################################################################################################################

- Script to plot the source radial profile+PSF profile and the residual profile for all wavelengths for each source.

- Plotting results derived from the script described above. 

###### Input files required ####################
1) .csv table files containing source radial profile, psf profile and residual profile data for each wavelength. eg:For CIT6 - cit6_UnInterpolated_70.csv, cit6_UnInterpolated_160.csv, cit6_UnInterpolated_450.csv, cit6_UnInterpolated_850.csv.
2) .csv table containing source distances in order to convert from projected to physical radius. - Source_Distances.csv


############ Output ##############

- plot of source radial profile+PSF profile and the residual profile for all wavelengths for the chosen source. Radial (x) axes are given in both projected (arcsec) and physical (pc) sizes.

"""

font = {'family' : 'normal',
        'size'   : 18,
	'style'  : 'oblique',
	'weight' : 'medium'}

plt.rc('font', **font)

Source_Data = Table.read('Source_Distances.csv', format='ascii.csv')

for Source in Source_Data:

	star = Source['Source']
	distance = Source['Distance(pc)']  #Units = pc


	#Distance for converion from arcsec to pc in the plots 

	#Loading and opening PACS 70 micron data files
	Stellar_Data_70 = Table.read(star+'_UnInterpolated_70.csv', format='ascii.csv')
	x_70 = Stellar_Data_70['x_axis(arcsec)']
	#x_70_cm = (distance * x_70) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	x_70_pc = (distance * x_70) * 4.8481e-6 #(pc * arcsec) is in AU's. So * 4.8... to convert to pc #4.8481e-6
	stellar_70 = Stellar_Data_70['Stellar_Profile(Jy/arcsec^2)']
	stellar_70 = stellar_70 * 1000 #Converting from Jy to mJy
	stellar_unc_70 = Stellar_Data_70['Stellar_Unc(Jy/arcsec^2)']
	stellar_unc_70 = stellar_unc_70 * 1000 #Converting from Jy to mJy
	psf_70 = Stellar_Data_70['PSF_Profile(Jy/arcsec^2)']
	psf_70 = psf_70 * 1000 #Converting from Jy to mJy
	res_70 = Stellar_Data_70['Residual_Profile(Jy/arcsec^2)']
	res_70 = res_70 * 1000 #Converting from Jy to mJy
	res_unc_70 = Stellar_Data_70['Residual_Unc(Jy/arcsec^2)']
	res_unc_70 = res_unc_70 * 1000 #Converting from Jy to mJy

	#Loading and opening PACS 160 micron data files
	Stellar_Data_160 = Table.read(star+'_UnInterpolated_160.csv', format='ascii.csv')
	x_160 = Stellar_Data_160['x_axis(arcsec)']
	x_160_cm = (distance * x_160) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	stellar_160 = Stellar_Data_160['Stellar_Profile(Jy/arcsec^2)']
	stellar_160 = stellar_160 * 1000 #Converting from Jy to mJy
	stellar_unc_160 = Stellar_Data_160['Stellar_Unc(Jy/arcsec^2)']
	stellar_unc_160 = stellar_unc_160 * 1000 #Converting from Jy to mJy
	psf_160 = Stellar_Data_160['PSF_Profile(Jy/arcsec^2)']
	psf_160 = psf_160 * 1000 #Converting from Jy to mJy
	res_160 = Stellar_Data_160['Residual_Profile(Jy/arcsec^2)']
	res_160 = res_160 * 1000 #Converting from Jy to mJy
	res_unc_160 = Stellar_Data_160['Residual_Unc(Jy/arcsec^2)']
	res_unc_160 = res_unc_160 * 1000 #Converting from Jy to mJy

	#Loading and opening SCUBA2 450 micron data files
	Stellar_Data_450 = Table.read(star+'_UnInterpolated_450.csv', format='ascii.csv')
	x_450 = Stellar_Data_450['x_axis(arcsec)']
	x_450_cm = (distance * x_450) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	stellar_450 = Stellar_Data_450['Stellar_Profile(mJy/arcsec^2)']
	stellar_unc_450 = Stellar_Data_450['Stellar_Unc(mJy/arcsec^2)']
	psf_450 = Stellar_Data_450['PSF_Profile(mJy/arcsec^2)']
	res_450 = Stellar_Data_450['Residual_Profile(mJy/arcsec^2)']
	res_unc_450 = Stellar_Data_450['Residual_Unc(mJy/arcsec^2)']

	#Loading and opening SCUBA2 850 micron data files
	Stellar_Data_850 = Table.read(star+'_UnInterpolated_850.csv', format='ascii.csv')
	x_850 = Stellar_Data_850['x_axis(arcsec)']
	x_850_cm = (distance * x_850) * 1.496e+13 #(pc * arcsec) is in AU's. So * 1.4... to convert to cm
	stellar_850 = Stellar_Data_850['Stellar_Profile(mJy/arcsec^2)']
	stellar_unc_850 = Stellar_Data_850['Stellar_Unc(mJy/arcsec^2)']
	psf_850 = Stellar_Data_850['PSF_Profile(mJy/arcsec^2)']
	res_850 = Stellar_Data_850['Residual_Profile(mJy/arcsec^2)']
	res_unc_850 = Stellar_Data_850['Residual_Unc(mJy/arcsec^2)']



	#################### Plotting ##################################################################


	fig = plt.figure(figsize=(16, 13)) #(width,height)
	gs = gridspec.GridSpec(4, 2) #(no.of.row, no.of.cols, height.ratio.of.rows)
	gs.update(hspace=0)
	gs.update(wspace=0)

	#Adding common x and y labels
	fig.text(0.52, 0.055,"Radius ($^{\prime\prime}$)", ha='center') #(distance.from.left.edge.of.image, distance.from.bottom.of.image)
	#fig.text(0.52, 0.93, 'Radius ($\\times 10^{18}$cm)', ha='center') #cm units
	fig.text(0.52, 0.93, 'Radius (pc)', ha='center') #pc units
	fig.text(0.06, 0.5,"log(Surface Brightness (mJy/arcsec$^2$))", va='center', rotation='vertical')
	fig.text(0.95, 0.5,"Residual (mJy/arcsec$^2$)", va='center', rotation=-90)


	#Plotting 70micron Radial + PSF Profile
	ax1 = fig.add_subplot(gs[0,0])
	radial_limit_70 = x_70 < 140
	ax1.errorbar(x_70[radial_limit_70], stellar_70[radial_limit_70], yerr=stellar_unc_70[radial_limit_70], fmt='o', color='midnightblue', markersize=3, label='70$\\mu m$ Stellar Profile') 
	ax1.errorbar(x_70[radial_limit_70], psf_70[radial_limit_70], fmt='-', markersize=0.5, mew=5, color='grey', linewidth=2, label='70$\\mu m$ PSF', alpha=0.6) 
	nbins = len(ax1.get_yticklabels())
	ax1.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
	ax1Ticks = ax1.get_xticks() 
	plt.setp(ax1.get_xticklabels(), visible=False)
	ax1.set_yscale('log')
	ax1.set_ylim([10**(np.log10(max(stellar_70[radial_limit_70])) - 5), 10**(np.log10(max(stellar_70[radial_limit_70])))])
	ax1.legend(fontsize=11)


	#Plotting 70micron Residual Profile
	ax2 = fig.add_subplot(gs[0,1], sharex=ax1) 
	radial_limit_70[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res profile
	ax2.errorbar(x_70[radial_limit_70], res_70[radial_limit_70], yerr=res_unc_70[radial_limit_70], fmt='^', color='orangered', markersize=3, mec='orangered', label='70$\\mu m$ Residual')
	ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
	ax2.yaxis.tick_right()
	ax2.yaxis.set_label_position("right")
	plt.setp(ax2.get_xticklabels(), visible=False)
	nbins = len(ax2.get_yticklabels())
	ax2.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max and min val on y axis to prevent merging with upper labels
	ax2.set_ylim(ymin=0) #starting y axis from zero
	ax2.legend(fontsize=11, loc=1)


	#Plotting 160micron Radial + PSF Profile
	ax3 = fig.add_subplot(gs[1,0], sharex=ax1)
	radial_limit_160 = x_160 < 130
	ax3.errorbar(x_160[radial_limit_160], stellar_160[radial_limit_160], yerr=stellar_unc_160[radial_limit_160], fmt='o', color='midnightblue', markersize=3,  label='160$\\mu m$ Stellar Profile') 
	ax3.errorbar(x_160[radial_limit_160], psf_160[radial_limit_160], fmt='-', markersize=0.5, mew=5, color='grey', linewidth=2, label='160$\\mu m$ PSF', alpha=0.6) 
	nbins = len(ax3.get_yticklabels())
	ax3.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
	plt.setp(ax3.get_xticklabels(), visible=False)
	ax3.set_ylim([10**(np.log10(max(stellar_160[radial_limit_160])) - 4.5), 10**(np.log10(max(stellar_160[radial_limit_160])))])
	ax3.set_yscale('log')
	ax3.legend(fontsize=11, loc=1)


	#Plotting 160micron Residual Profile
	ax4 = fig.add_subplot(gs[1,1], sharex=ax1) 
	radial_limit_160[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res profile
	ax4.errorbar(x_160[radial_limit_160], res_160[radial_limit_160], yerr=res_unc_160[radial_limit_160], fmt='^', color='orangered', mec='orangered',  markersize=3, label='160$\\mu m$ Residual')
	ax4.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
	ax4.yaxis.tick_right()
	ax4.yaxis.set_label_position("right")
	plt.setp(ax4.get_xticklabels(), visible=False)
	nbins = len(ax4.get_yticklabels())
	ax4.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
	ax4.set_ylim(ymin=0) #starting y axis from zero
	ax4.legend(fontsize=11, loc=1)


	#Plotting 450micron Radial + PSF Profile
	ax5 = fig.add_subplot(gs[2,0], sharex=ax1)
	radial_limit_450 = x_450 < 70
	ax5.errorbar(x_450[radial_limit_450], stellar_450[radial_limit_450], yerr=stellar_unc_450[radial_limit_450], fmt='o', color='midnightblue', markersize=3,  label='450$\\mu m$ Stellar Profile') 
	ax5.errorbar(x_450[radial_limit_450], psf_450[radial_limit_450], fmt='-', markersize=0.5, mew=5, color='grey', linewidth=2, label='450$\\mu m$ PSF', alpha=0.6) 
	nbins = len(ax5.get_yticklabels())
	ax5.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both'))
	plt.setp(ax5.get_xticklabels(), visible=False)
	ax5.set_ylim([10**(np.log10(max(stellar_450[radial_limit_450])) - 5), 10**(np.log10(max(stellar_450[radial_limit_450])))])
	ax5.set_yscale('log') #setting y axis to log scale
	ax5.legend(fontsize=11, loc=1)


	#Plotting 450micron Residual Profile
	ax6 = fig.add_subplot(gs[2,1], sharex=ax1) 
	radial_limit_450[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res profile
	ax6.errorbar(x_450[radial_limit_450], res_450[radial_limit_450], yerr=res_unc_450[radial_limit_450], fmt='^', color='orangered', mec='orangered', markersize=3, label='450$\\mu m$ Residual')
	ax6.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
	ax6.yaxis.tick_right()
	ax6.yaxis.set_label_position("right")
	plt.setp(ax6.get_xticklabels(), visible=False)
	nbins = len(ax6.get_yticklabels())
	ax6.yaxis.set_major_locator(MaxNLocator(prune='both'))
	ax6.set_ylim(ymin=0) #starting y axis from zero
	ax6.legend(fontsize=11, loc=1)


	#Plotting 850micron Radial + PSF Profile
	ax7 = fig.add_subplot(gs[3,0], sharex=ax1)
	radial_limit_850 = x_850 < 90
	ax7.errorbar(x_850[radial_limit_850], stellar_850[radial_limit_850], yerr=stellar_unc_850[radial_limit_850], fmt='o', color='midnightblue', markersize=3, label='850$\\mu m$ Stellar Profile') 
	ax7.errorbar(x_850[radial_limit_850], psf_850[radial_limit_850], fmt='-', markersize=0.5, mew=5, color='grey', linewidth=2, label='850$\\mu m$ PSF', alpha=0.6) 
	nbins = len(ax7.get_yticklabels())
	ax7.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper'))
	ax7.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels
	ax7.set_ylim([10**(np.log10(max(stellar_850[radial_limit_850])) - 5), 10**(np.log10(max(stellar_850[radial_limit_850])))])
	ax7.set_yscale('log') #setting y axis to log scale
	ax7.legend(fontsize=11, loc=1)


	#Plotting 850micron Residual Profile
	ax8 = fig.add_subplot(gs[3,1], sharex=ax1) 
	radial_limit_850[0:2]=False #We make the 0th and 1st point a false index since we don't want to plot the 1st two points in the res profile
	ax8.errorbar(x_850[radial_limit_850], res_850[radial_limit_850], yerr=res_unc_850[radial_limit_850], fmt='^', color='orangered', mec='orangered', markersize=3, label='850$\\mu m$ Residual')
	ax8.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
	ax8.yaxis.tick_right()
	ax8.yaxis.set_label_position("right")
	nbins = len(ax8.get_xticklabels())
	ax8.yaxis.set_major_locator(MaxNLocator(prune='both')) #removing max val on y axis to prevent merging with upper axis labels
	ax8.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels
	ax8.set_ylim(ymin=0) #starting y axis from zero
	ax8.legend(fontsize=11, loc=1)


	#Adding second x axis on top of plot for the cm x axis scale using the cm scale for 70micron for both stellar and res profiles.
	ax9 = ax1.twiny()
	#x_70_cm = x_70_cm / 1e18
	#ax9.errorbar(x_70_cm, stellar_70, yerr=stellar_unc_70, fmt='o', color='none', mec='none') #cm units
	ax9.errorbar(x_70_pc, stellar_70, yerr=stellar_unc_70, fmt='o', color='none', mec='none') #pc units
	ax9.set_yscale('log') #setting y axis to log scale
	ax9.set_ylim([10**(np.log10(max(stellar_70)) - 5), 10**(np.log10(max(stellar_70)))])
	ax9.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels


	ax10 = ax2.twiny()
	#ax10.errorbar(x_70_cm, res_70, yerr=res_unc_70, fmt='^', color='none', mec='none') #cm units!
	ax10.errorbar(x_70_pc, res_70, yerr=res_unc_70, fmt='^', color='none', mec='none') #pc units!
	ax10.xaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='both')) #removing max val on x axis to prevent merging with res axis labels
	ax10.set_ylim(ymin=0) #starting y axis from zero

	save_file = star + '_RadialProfile.png'
	plt.savefig(save_file, bbox_inches='tight')
	plt.show()




















