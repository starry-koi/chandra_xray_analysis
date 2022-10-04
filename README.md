Chandra X-ray Data Analysis Pipeline
===

Analysis of Chandra X-ray images via python code with CIAO integration. Finds the X-ray flux (in a specific band) for sources within/around a galaxy, following along with the procedure regarding X-ray data in [this 2021 paper](https://arxiv.org/abs/2105.05876). Currently only works if the galaxy you want to analyze is located at the aimpoint of the S3 chip (i.e. if it was the target of the observation). The coordinates of the aimpoint can be checked in the pdf included with the observation, assuming the data was downloaded from the Chandra archive.

Be sure to read the 'Rerunning the Code' section - there's an important warning there that if not heeded could mess up your analysis.

Also -- currently there is an error with CIAO's `wcs_match` function (used in the code to fix the astrometry) that causes a segfault when there are no matches between the source list and catalog source list. To (kind of) get by this, simply add the obsid of the problem galaxy to the `skip_obsid_astrom` variable *as a string* and it should skip right over the use of the `wcs_match` function for that galaxy. Note that this means the astrometry won't be changed for the galaxy, though since the bug should only be happening when there were no matches between our found sources and the catalog sources, this should be alright.


Setup
---

You will need to download and install the current version of [CIAO](https://cxc.cfa.harvard.edu/ciao/) (see the detailed installation instructions on their site). This script currently supports CIAO 4.14, which is the latest version as of time of writing. You will also need to place the folder containing `xray_flux.py` in the same directory as the individual OBSID folders for the observations you wish to analyze.


Brief Rundown of the Code
---

A more in-depth walkthrough is outlined in the header in the `xray_flux.py` code file, and details about specific sections and steps are detailed in the comments throughout the code. Here we provide a brief overview.

* We assume that the galaxy of interest is the target of the Chandra observation. We find the evt2 file in the OBSID folder and reprocess the data using the CIAO `chandra_repro` function, then restrict further analysis to the S3 chip (as that is the chip the aimpoint is located on). This creates a /repro folder in the OBSID folder with the labelling schema /repro_\[band_check]-band__\[filter_check]-filter; e.g. if you're filtering from 2-7 keV and finding flux in the 2-10 keV band, the folder name will be /repro_2-10-band__2-7-filter. Any further created files and folders will be located in this /repro folder, save for the final results folder which is located in the code folder (with `xray_flux.py`). 

* We correct the astrometry of the Chandra image by running the CIAO functions `fluximage` and `wavdetect` on the Chandra image to obtain a list of sources and comparing these to a retrieved list of SDSS catalog sources. 

* We next filter the corrected image for background flares.

* We start our analysis on the cleaned, corrected Chandra image by filtering it to whatever band is specified by the `filter_check` variable and running `fluximage` and `wavdetect` to obtain a list of detected sources. We then restrict our interest to sources falling inside the region of the galaxy we care about, and find the source fluxes in the energy band specified by the `band_check` variable using the CIAO function `srcflux`.

* We calculate the errors on the detected source net counts and positions, etc.

* Finally, we save the results in a created results folder specific to the filtering band and flux band (`filter_check` and `band_check` variables, respectively). This folder is located in the code folder.


Variables the user will have to care about
---

You will have to go into the code file to change the values of at least some of these variables depending on what galaxies you are analyzing and how you are analyzing them. Comments in `xray_flux.py` give more detail, with an overview given here.


* `r50_all_gals` 
	* The Petrosian 50% light radius for each galaxy, in order of ascending OBSID. This, along with the `gal_exclude_rad` variable, is used to define the 'area of interest' around each galaxy that we'll use to decide what sources we want to analyze.

* `gal_exclude_rad` 
	* The galaxy region (that 'area of interest') is defined by taking each galaxy's specific r50 value and multiplying it by this. Currently it's set to 3, such that the galaxy region we focus on has a radius of 3\*r50, which seems to work well.

* `galdist` and `galdist_flag` 
	* Used to automatically find luminosities from calculated source fluxes, the first is a list of the distances to each galaxy (in Mpc) while the second tells the code whether or not you want it to find luminosities.

* `band_check`
	* Which energy band (in keV) you want to find source fluxes in; e.g. 0.5-2, 2-10, etc.

* `filter_check`
	* Which energy band (in keV) you want to filter the cleaned, corrected image by before finding and analyzing X-ray sources.


There are a number of other variables detailed in the code that can be changed by the user, but I've found that they rarely need to be changed from galaxy to galaxy and project to project.


A Sample Walkthrough of a Code Run
---

Download your OBSIDs from the [Chandra Data Archive](https://cda.harvard.edu/chaser/) and make sure that the galaxies you want to analyze were the actual targets of the Chandra observations by checking the RA and Dec in the pdfs in each OBSID folder. Place the folder containing `xray_flux.py` in the same directory as the individual OBSID folders. Go into the `xray_flux.py` script and update the `r50_all_gals` and `galdist` variables, making sure the `galdist_flag` variable is set to True. Check that the `band_check` and `filter_check` variables are set to the desired energy ranges.

Once everything is set up, open a terminal in the code directory and initialize CIAO, then run the code as

	$./xray_flux.py

As the code runs, it will print out where it is in the analysis process so you know it's working and what it's doing.


Understanding the Output
---

While the code makes and saves a number of intermediary files during analysis, these are all saved under the various /repro directories in each OBSID folder (see the 'Rerunning the Code' section for more detail). Everything that you're likely to care about is saved to the /results folder in the code directory - this folder is labelled according to the `band_check` and `filter_check` variables, similar to the /repro folders. Explanations of what each file contains is as follows:
* `galaxy_data.txt`
	* information that only changes from galaxy to galaxy - e.g. astrometry-related information, how many good/bad sources there were, etc.
* `params.txt`
	* the specific values of the parameters used in the run of the code. See the warning in this file (and under the 'Rerunning the Code' section in this ReadMe) about changing time-intensive-task related parameters without deleting the corresponding folders.
* `src_all.txt`
	* all the information regarding detected sources - counts, net counts, information on source and background apertures, etc.
* `src_summary.txt`
	* handpicked summary of the important information - net counts, fluxes, etc.
* /cleaned_lightcurves
	* images of the cleaned lightcurves for each galaxy (after removing background flares)
* /regions_and_images
	* Each galaxy will have a number of files associated with it (some depending on if any X-ray sources were detected in the galaxy region) as detailed below:
	* `fin_psf_srcs.bg.reg`
		* the background regions used for the sources. Note that when a background region intersected with a source region, the intersection was excluded from the background region in calculations. (This region will only display correctly when used on the Chandra .fits file associated with the observation, such as the `new_evt2.fits` file).
	* `fin_psf_srcs.src.reg`
		* the source regions used to calculate the fluxes. Note that when multiple source regions intersected, the intersection was excluded from both source regions for calculations. (This region has the same display requirements as `fin_psf_srcs.bg.reg`).
	* `galaxy_region_ds9.reg`
		* the region used to define 'in the galaxy' - fluxes are only calculated for sources inside this region.
	* `new_evt2.fits`
		* the filtered, cleaned, and astrometrically corrected X-ray images used to calculate fluxes
	* `xray_ds9_regions.reg`
		* X-ray source positions with positional errors and names (in RA, Dec).


Rerunning the Code
---

If you are rerunning the code on the same galaxies with the same `band_check` and `filter_check` variables (but with other variables changed) you will need to delete the relevant /repro folder in the OBSID folder(s) and run the code fresh. This is because for various time-intensive tasks (`chandra_repro`, `fluximage`, `wavdetect`, and `srcflux`) the code checks to see if they've already been run by looking for specific folders. If it finds the folders it assumes the task has already been run and skips it. This means that if you change one of the variables related to one of those tasks, and don't delete the folder that the code is checking for, it will assume that the task has already been run and use the previous results (without using the updated variable) *and* it won't display a warning or anything. So this is your warning. 

Folder names and what task they correspond to (the latter three are located in the /repro folder):
* /repro_...
	* `chandra_repro` - deleting this folder will have the code restart completely fresh.
* /fixastrom_nogal
	* astrometry-related `fluximage` and `wavdetect` - deleting this will force the code to redo the astrometry.
* /minflux_calc
	* `srcflux` used for calculating minimum fluxes for galaxies.
* /analysis
	* analysis-related `fluximage` and `wavdetect` - used to find sources for analysis - deleting this will force the code to re-find the sources.
* /flux
	* `srcflux` used for finding the fluxes from the detected sources.

In general, best practice is that if you change any of the variables, best to just delete all the code-created folders (which can be done by deleting the relevant /repro folder(s) and the relevant /results folder) and run everything fresh. 


Quirks
---

Most of these are detailed elsewhere in the ReadMe, but here's a quick gathered list.

* Generally, the only thing you will need to change when using the code on different galaxies are the `r50_all_gals` and `galdist` variables.
* Similarly, the only thing you will need to change when using the code with a different filter or flux band are the `band_check` and `filter_check` variables.
* The code only works if the galaxy of interest is at the aimpoint of the S3 chip (i.e. if it was the target of the observation).
* The astrometry won't get updated correctly if there is more than one ASOL file attached to the .fits file. This would be a fairly quick fix in the code if it were to occur, but first I'd need to know the keywords for the extra ASOL files, which I don't. I haven't encountered any observations thus far that have multiple ASOL files, so this is unlikely to occur in the first place, though.
* The number of expected background sources are only calculated for hard and soft band flux bands (`band_check` set to 2-10 and 0.5-2, respectively), as those are the only bands the relevant equation is defined for. Other bands will have NaN in this column in the results folder.
* During the astrometry correcting, DS9 will pop up and open up a catalog. This is normal and supposed to happen. 
* Sometimes during this DS9 will freeze (at least on my machine) and I've found that pressing the alt key will get it to unfreeze. No idea why, though.
* The code is highly commented, so cruising through it might help you understand more what's going on.
* Run times per galaxy can take 3-10 minutes, depending on how fast your computer is and whether there are any X-ray sources found in the galaxy or not.


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
