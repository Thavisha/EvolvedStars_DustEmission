###################################################################################################################################
    ########### Radial Profile Generation ################################
###################################################################################################################################

- Generating radially avaraged surface brightness profiles using the SCUBA-2 and PACS observations. 

- The PSF profile of each instrument and each wavelength is used to derive the PSF-reduced residual profile. This will give the profile for the extended CSE compoenent.
- PACS PSFs are given in .fits file form (Bocchio et al., 2016) and the SCUBA-2 PSF is a Guassian with two componenets for the main beam and secondary beam (Dempsey et al., 2013)
- The central pixel (brightest pixel) were initally measured using CalTech Aperture Photometry Tool (APT), which fits a 2D Guassian to the source and finds its peak to find the peak position.
- Both interpolated (required for SED fitting) and uniterpolated (required for plotting) profiles are generated. 
 

###### Input ####################

1) .fits files of each wavelength for the chosen source. 
2) PACS PSF .fits files (one for each wavelength) downloaded from Boccio et al., 2016 via vizier
3) Source_Coordinate_List.csv file containing source central positions. 


###### Output ##################################  

1) PSF flux - as terminal print out
2) Characteristic extension (extension at peak point of stellar profile) - as terminal print out
3) Extension at maximum 3sigma brightness level - as terminal print out

4) .csv file for each source containing the radius, residual profile and +/-unc of residual profile  -- required for SED fitting
5) .csv file for each source contaiing the uniterpolated source profile and +/-unc of source profile, PSF profile, residual profile and +/-unc of residual profile, fractional residual profile and +/-unc of fractional residual profile -- required for plotting these profiles. (fractional residual profiles are not required for the paper. They are just extra information for testing).  


** Folder "Dharmawardena_RadialProfiles" give the radial profiles generated using this code and used in this paper. 



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





























