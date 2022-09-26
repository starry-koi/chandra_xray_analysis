Chandra X-ray Data Analysis Pipeline
===

Python code with CIAO integration for analysis of Chandra X-ray images. Finds the X-ray flux (in a specific band) for sources within/around a galaxy. Currently only works if the galaxy you want to analyze is located at the aimpoint of the S3 chip (i.e. if it was the target of the observation). The coordinates of the aimpoint can be checked in the pdf included with the observation (assuming the data was downloaded from the Chandra archive).


Setup
---

You will need to download and install the current version of [CIAO](https://cxc.cfa.harvard.edu/ciao/) (see the detailed installation instructions on their site). This script currently supports CIAO 4.14, which is the latest version as of time of writing.


Brief Rundown of the Code
---

A more in-depth walkthrough is outlined in the header in the `xray_flux.py` code file, and details about specific sections and steps are detailed in the comments throughout the code. Here we provide a brief overview.

We assume that the galaxy of interest is the target of the Chandra observation. We find the evt2 file in the OBSID folder and reprocess the data using the CIAO `chandra_repro` function, then restrict further analysis to the S3 chip (as that is the chip the aimpoint is located on). 

We correct the astrometry of the Chandra image by running the CIAO functions `fluximage` and `wavdetect` on the Chandra image to obtain a list of sources and comparing these to a retrieved list of SDSS catalog sources. 

We next filter the corrected image for background flares.

We start our analysis on the cleaned, corrected Chandra image by filtering it to whatever band is specified by the `filter_check` variable and running `fluximage` and `wavdetect` to obtain a list of detected sources. We then restrict our interest to sources falling inside the region of the galaxy we care about, and find the source fluxes in the energy band specified by the `band_check` variable using the CIAO function `srcflux`.

We calculate the errors on the detected source net counts and positions, etc.

Finally, we save the results in a created results folder specific to the filtering band and flux band (`filter_check` and `band_check` variables, respectively).


Variables the user will have to care about
---

You will have to go into the code file to change the values of at least some of these variables depending on what galaxies you are analyzing and how you are analyzing them. Comments in `xray_flux.py` give more detail, with an overview given here.


`r50_all_gals` - the Petrosian 50% light radius for each galaxy, in order of ascending OBSID. This, along with the `gal_exclude_rad` variable, is used to define the 'area of interest' around each galaxy that we'll use to decide what sources we want to analyze.

`gal_exclude_rad` - The galaxy region (that 'area of interest') is defined by taking each galaxy's specific r50 value and multiplying it by this. Currently it's set to 3, such that the galaxy region we focus on has a radius of 3\*r50, which seems to work well.

`galdist` and `galdist_flag` - Used to automatically find luminosities from calculated source fluxes, the first is a list of the distances to each galaxy (in Mpc) while the second tells the code whether or not you want it to find luminosities.

`band_check` - which energy band (in keV) you want to find source fluxes in; e.g. 0.5-2, 2-10, etc.

`filter_check` - which energy band (in keV) you want to filter the cleaned, corrected image by before finding and analyzing X-ray sources.


There are a number of other variables detailed in the code that can be changed by the user, but I've found that they rarely need to be changed from galaxy to galaxy and project to project.


A Sample Walkthrough
---

Download your OBSIDs from the [Chandra Data Archive](https://cda.harvard.edu/chaser/) and make sure that the galaxies you want to analyze were the actual targets of the Chandra observations by checking the RA and Dec in the pdfs in each OBSID folder. Place the folder containing `xray_flux.py` in the same directory as the individual OBSID folders. Go into the `xray_flux.py` script and update the `r50_all_gals` and `galdist` variables, making sure the `galdist_flag` variable is set to True. Check that the `band_check` and `filter_check` variables are set to the desired energy ranges.

Once everything is set up, open a terminal in the code directory and initialize CIAO, then run the code as

	$./xray_flux.py

As the code runs, it will print out where it is in the analysis process so you know it's working and what it's doing.

Understanding the Output
---


Rerunning the Code
---

The code puts everything OBSID-specific (other than final results) in a file folder in the OBSID folder labelled /repro_\[band_check]-band__\[filter_check]-filter; e.g. if you're filtering from 2-7 keV and finding flux in the 2-10 keV band, the folder name will be /repro_2-10-band__2-7-filter. Rerunning the code after changing the `band_check` and/or `filter_check` variables will create a new /repro folder with the new flux/filter bands.


In general, best practice is that if you change anything 

When 
In general, best practice is that if you change anything, best to delete the entire /repro folder in the OBSID folder corresponding to that run of the code.

The code checks whether time-intensive tasks (such as `chandra_repro`, `fluximage`, `wavdetect`, and `srcflux`) have been run by checking to see if the relevant folder exists. If the folder does exist, then the code won't run those functions again. So if you have already run the code and would like to change some variables regarding those tasks and have the code redo that analysis, 


Quirks
---


If something went wrong
---

* The code just stopped around the "Matching wavdetect and catalog sources" step
	* Currently there's a bug with the CIAO `wcs_match` function and sometimes it fails and throws out a segfault that stops the code. This should only be happening when there were no matches between our found sources and the catalog sources. I've implemented a workaround - you can put the OBSID of the problem galaxy to the `skip_obsid_astrom` variable (as a string) to tell the code to skip running the `wcs_match` function.
* I got an 'out of bounds' error for the `r50_all_gals` variable 
	* Likely forgot to update the r50 values to match with the OBSIDs you downloaded. Need to have one r50 value per OBSID.
* DS9 popped up and froze and the code timed out
	* Sometimes DS9 just freezes, which is annoying. I've found that pressing the alt key gets DS9 to unfreeze so the code can proceed.
* Error popped up saying something was wrong with 'from ciao_contrib.runtool import \*'
	* Likely forgot to initialize CIAO before running the code. 
* There are NaN values in the "expected background sources" column of the text files in the results folder
	* The expected background sources can only be calculated for soft (0.5-2 keV) and hard (2-10 keV) bands, so if you're using a different band then this column will have NaNs in it as no value could be calculated.
* The `fin_psf_srcs.bg.reg` and `fin_psf_srcs.src.reg` files aren't displaying correctly in DS9
	* Those region files are in pixel coordinates, so they will only display correctly when used on the Chandra image associated with that observation (such as the `new_evt2.fits` file).
* The `fin_psf_srcs.bg.reg` and `fin_psf_srcs.src.reg` files have regions that look like they intersect, is the flux from the intersecting part being used in both sources/backgrounds?
	* No, overlapping parts of source/background region files are excluded from calculations.





# xray_code
Python code with CIAO integration for analysis of Chandra X-ray images

--- Currently there is an error with CIAO's wcs_match function (used in the code to fix the astrometry) that causes a segfault when there are no matches between the source list and catalog source list. To (kind of) get by this, simply add the obsid of the problem galaxy to the skip_obsid_astrom variable AS A STRING and it should skip right over the use of the wcs_match function.


GENERAL INFO--------------
In order for this code to work you'll need to manually enter the Petrosian 50% light radius (r50) value for each galaxy - see the VARIABLES section in the code. You also have the option to add in the galaxy distances to calculate luminosities - if you do be sure to set galdist_flag to True to calculate them (see VARIABLES section). You will also need to specify the energy bands for filtering the final X-ray image and the energy bands to be used for calculating flux (see the VARIABLES section).

You should put the folder containing this code in the directory containing all the obsid folders for your observations. Then, after changing the values of r50_all_gals (and galdist and galdist_flag if you want to calculate luminosities), as well as band_check and filter_check, you should just be able to run the code and it'll work. Note that the code checks to see whether time-intensive tasks (such as chandra_repro, fluximage, wavdetect, and srcflux) have already been run by checking to see if the folder exists - and if the folder exists it won't run it again. So if you have already run the code, and want to change parameters regarding those functions and have the code redo the functions, you'll have to delete the relevant folder (in the obsid directory):
-repro_.../		-  chandra_repro  (note doing this will basically have the code restart completely from zero)
-fixastrom_nogal/	-  astrometry-related fluximage and wavdetect (used to correct the astrometry)
-analysis/		-  analysis-related fluximage and wavdetect (used to find sources for analysis) 
-flux/			-  srcflux (finding the actual fluxes from the sources)

Note that the repro/ folder will also have the flux band and filtering band in the title - e.g., if you're filtering the image at 2-7 keV and finding flux from 2-10 keV, the folder name will be "repro_2-10-band__2-7-filter". Also important - the code creates a results folder to store results in - the name varies depending on the specific band and filter, in exactly the same way as the repro folder. Before rerunning the code, you should also delete the relevant results folder generated by the first run of the code, or else things might end up messy with duplicate files and such. This can also be avoided by renaming that initial results folder e.g. to /results_old


QUIRKS/TIPS---------------
-The only thing you need to change when using the code on different galaxies is the r50_all_gals and galdist variables.
-The only thing you need to change when using the code with a different filter or flux band is the band_check and filter_check variables.
-The code only works if the galaxy of interest is at the aimpoint of the S3 chip (i.e. if it was the target of the observation).
-The astrometry won't get updated correctly if there is more than one ASOL file. This would be a fairly quick fix in the code if it were to occur, though.
-The expected background sources are only calculated for hard and soft band fluxes (2-10 and 0.5-2 keV, respectively), as those are the only bands the relevant equation is defined for. Other bands will have NaN in this column in the results folder.
-Sometimes DS9 freezes (at least on my machine) e.g. while doing the astrometry matching. You can press the alt key to get it to unfreeze. No idea why.
-Read the entire readme file and at least the overview in the script, as they answer most questions. The code is also highly commented, so cruising through that might help you solve your issue as well.


OUTPUT--------------------
Everything of importance will go in the results folder. 
-galaxy_data.txt	-  information that only changes from galaxy to galaxy - e.g. astrometry-related information, how many good/bad sources there were
-params.txt		-  the specific values of the parameters used in the run of the code. SEE THE WARNING IN THIS FILE about changing time-intensive-task related parameters without deleting the corresponding folders.
-src_all.txt		-  all the information regarding specific sources - counts, net counts, information on source and backround apertures, etc.
-src_summary.txt	-  handpicked summary of the important information - net counts, fluxes, etc.

The regions_and_images folder  -  This one needs some extra explanation. Multiple files are included for each galaxy that had sources found in it:
-fin_psf_srcs.bg.reg	-  the background regions used for the sources. Note that when a background region intersected with a source region, the intersection was excluded from the background region in calculations.
-fin_psf_srcs.src.reg	-  the source regions used to calculate the fluxes. Note that when multiple source regions intersected, the intersection was excluded from both source regions for calculations.
	A note about these two region files - they will only work when used on the Chandra .fits file associated with that observation (such as the new_evt2.fits file) as the region file coords are in pixels.
-galaxy_region_ds9.reg	-  the galaxy region used - sources inside this region were excluded when matching found sources to catalog sources, and we only calculated fluxes for sources inside this region.
-new_evt2.fits		-  the X-ray image (filtered, cleaned, and astrometrically corrected) used to calculate fluxes
-xray_ds9_regions.reg	-  source positions with positional errors and names (in sky coordinates - use these if you want to put the sources on e.g. an HST image)



OTHER--------------------
The code follows along with the procedure regarding X-ray data in [this 2019 paper](https://arxiv.org/abs/1907.12585) and [this 2021 paper](https://arxiv.org/abs/2105.05876). We use Dickey & Lockman maps for our galactic column densities, and we use a power-law spectral model with a photon index of 1.8. For a more in-depth explanation, see the code. There's a brief overview at the top, and the code is divided into sections - under each section is a slightly more in-depth explanation of what's happening, and then the individual comments should explain in detail what's going on.

Common issues/quirks are explained in the overview in the code.

Rerunning - if you are rerunning the code on the same galaxies with the same flux and filtering bands (but with other parameters changed) you will need to delete the relevant repro/ folder and rerun the time-intensive tasks in order to get those parameters to actually take. In general, if you change any of the parameters, best practice is to delete the repro/ folder just in case and rerun the code from scratch.

The source and background final regions (fin_psf_srcs.bg.reg and fin_psf_srcs.src.reg) will only display correctly on the Chandra image of the galaxy (new_evt2.fits), as they're in pixel coordinates and not sky coordinates. e.g. if you try to overlay those regions on an HST image in DS9, the regions won't be in the correct place.


ADDING NEW BANDS--------------------

If you want to add a new option to band_check, you'll need to ctrl+f and add an additional elif statement where it's called in the code. Additionally, you'll need to add another option to srcflux_band, which means you'll need to find the appropriate specific energy to use for the new band, which can most likely be found on the CIAO website.


