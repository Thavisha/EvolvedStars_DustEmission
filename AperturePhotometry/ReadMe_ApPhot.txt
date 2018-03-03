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







































