#!/usr/bin/env python
from ciao_contrib.runtool import *
import numpy as np
import glob
import os
import regions as reg
import region #yes, this is different than the regions package. I know.
import poissonstats as ps #user-made file, has functions for finding the error for counts
import csv
import sys
import re
import math
import collections
import warnings
import pandas as pd
import subprocess as sp
from sh import gunzip
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from decimal import Decimal

#Author: Lilikoi Latimer
#Created October 2019
#Updated September 2022
#Currently works under CIAO 4.14

'''
======================================================================================================


### Quirks and possible issues ###

- Read the readme file as well as this header.

- Need to have the folder containing this code (folder xray_flux_code) in the same directory as the 
  obsids. Because what this code does is start in the xray_flux_code directory, then go up one
  directory, then go into a specific obsid directory.

- Only works if the galaxy we care about is at the aimpoint of the S3 chip i.e. the telescope was
  pointed right at the galaxy.

- Will need to MANUALLY set the bands and filtering bands (band_check and filter_check in the 
  VARIABLES section) to appropriate bands. filter_check will be used to filter the final image, and 
  band_check will be used to determine what energy band to calculate flux for.

- Have to MANUALLY set the Petrosian 50% light radii for each galaxy (r50_all_gals, see VARIABLES 
  section). Can also set the distance to the galaxies, in Mpc (galdist) and galdist_flag to True if
  you want to calculate the log luminosities as well. Everything else should run reasonable well if
  just left as-is (at least for finding AGN).

- Different band_check and filter_check will store files in different directories in the OBSID folders,
  appropriately labelled, so you will only need to delete folders and rerun code if you're trying to
  redo everything at the same band and filter as you previously did (see below).
  
- The code will skip doing time-intensive tasks (reprocessing, running fluximage, wavdetect, and
  srcflux) if it's detected that it's already done those before - it does this by searching for 
  directories that contain the outputs of those tasks, and only runs them again if it can't find the
  directory. --So if you want to run, say, srcflux again because you changed some parameters, then 
  you'll need to go into the relevant obsid/repro/ directory and delete the flux folder.
  To run fluximage,wavdetect for ASTROMETRY again -> delete fixastrom_nogal folder in obsid/repro/
  To run fluximage,wavdetect for ANALYSIS again   -> delete analysis folder in obsid/repro/
  To run srcflux again				  -> delete flux folder in obsid/repro/

- Related to the above, if you run the code with a set of parameters, then change the parameters for
  one of the resource intensive functions (fluximage,wavdetect,srcflux), then run the code again, the
  parameter output file will change (indicating that it ran with the new parameters) but since the
  code will skip actually doing those tasks, it won't have actually run the functions with the new
  parameters. See results/params.txt for an explanation, or scroll down to the bottom of the script
  to see it written out there.

- Code might bug out if the initial (before reprocessing) evt2 files (found under obsid/primary) are 
  zipped. If it does, just run the code again and the problem should go away. Or, alternatively, unzip 
  them manually.

- This code isn't written to handle a fits evt2 file that has multiple ASOL files attached. It would 
  be easy to change, but I'd need to know the keywords for multiple ASOL files in the header, and I
  don't, so it's just left as is currently.

- At the moment, wcs_match (the CIAO function) appears to be throwing errors and segfaults when there
  are zero matches between source list and catalog list. This can occasionally stop the code in the
  astrometry-fixing section. As a quick fix, I've thrown in a variable skip_obsid_astrom, which, if
  you put the obsid of the problem-causing galaxy in AS A STRING then should skip right over the 
  problem section and allow you to get fluxes and everything. Hopefully CIAO fixes this soon though.

### Broad Overview ###

This program finds the xray flux in a specific band for sources within/around a galaxy (the exact region depends on user input, see VARIABLES section). It works for multiple observations on multiple galaxies, as well.

First, we find the evt2 file in obsid/primary/ (unzipping the file if necessary), and from the header of the file we retrieve the galaxy's name, ra, and dec (note that we assume the galaxy is the target of the observation). We then reprocess the data, taking into account whether the observations were taken in FAINT or VFAINT mode. From here on out, we focus on the S3 chip, as we're assuming the galaxy is the target of the observation, which means that's the only chip we care about.

We fix the astrometry by running fluximage and wavdetect on the reprocessed evt2 file to get a list of sources (keeping only ones with significance greater than 5), which we then compare to a list of catalog sources (taken from the SDSS catalog). Note that when fixing the astrometry we exclude any wavdetect sources that lie within the galaxy region. We find the matches, and if the number of matches is greater than the minimum required (min_matches - see VARIABLES section) we go ahead and fix the astrometry, updating the fits file and the ASOL file.

Next, we take our astrometrically-corrected image and filter it for background flares (while doing so we exclude the galaxy region and any wavdetect sources we found prior - see FILTERING and CLEANING section).

We then start our analysis. We filter the corrected, cleaned image to whatever band is specified in filter_check (see VARIABLES section) and run fluximage and wavdetect. Afterwards, we focus in on only the wavdetect regions that fall within our galaxy region. We use the positions of these wavdetect regions as a basis for our source and background regions, which we use to find the flux via srcflux at the energy band specified by band_check (see VARIABLES section). We extract the data we're interested in from the srcflux output file, such as source positions, counts, background counts, etc.

Now we calculate our errors. We find errors on the net counts and source position, check to see whether the sources are significant or not, and find the expected number of background sources to fall within our galaxy region. This can get a bit complicated - see CALCULATING ERRORS section for more detail.

Finally, we save the results. We write a summary file on all the sources we've found, including only some of the data we pulled, like source positions, net counts, fluxes, (the important, interesting stuff) and also write out a file with all the data we calculated or pulled from srcflux. We also write comments in the files explaining what each column means. We write out a file concerning galaxy-specific information (was the astrometry corrected, by how much, what were the expected background counts in the galaxy region, etc.), and also a folder containing the final cleaned images and any relevant regions for DS9 (galaxy region, source regions) as well as a folder containing the cleaned lightcurves for our galaxies. And lastly, we write out a parameter file that shows the values of all the parameters used for the specific run.




======================================================================================================
'''


#=====================================================================================================
#=====================================================================================================

#############################################################
################# VARIABLES and DIRECTORIES #################
#############################################################
'''
Defaults of:

r50_all_gals 		= [specific to the galaxies]
galdist 		= [specific to the galaxies]
galdist_flag 		= True 
band_check		= [user input]
filter_check		= [user input]
coord_frame 		= 'icrs'
gal_exclude_rad 	= 3
fluxim_astrom_psfecf 	= 0.9
wavdet_astrom_scales 	= '1.0 1.4 2.0 2.8 4.0'
catalog_filter 		= '$imag<22'
astrom_match_radius 	= 3
astrom_match_residlim 	= 2
astrom_match_method 	= 'trans'
min_matches 		= 1
bin_time_var 		= 200
fluximage_spec_energy 	= 4
fluximage_psfecf 	= 0.393
wavdet_analysis_scales 	= '1.0 1.4 2.0 2.8 4.0'
analysis_psfecf 	= 0.9
analysis_energy 	= 4.5
bg_radius_factor 	= 12 
srcflux_PhoIndex	= '1.8' 
minflux_counts 		= 2
minflux_bgcounts 	= 0
ds9_reg_color 		= 'magenta'
round_dec 		= 2 
flux_ref 		= 1e15
mult_hdrs 		= False
'''


#all values ordered from lowest obsid to highest
#moved these variables up here from the galaxy-specific variables section for ease of use
r50_all_gals = [8.05, 8.75, 7.85]
galdist = [80.21, 122.78, 107.6]
galdist_flag = True

skip_obsid_astrom = [] #Due to wcs_match throwing errors, if any specific obsid makes the program
#stop due to a segfault during the astrometry fixing phase (specifically when matching sources
#to catalog sources) just put the obsid in this list AS A STRING and it should skip the step
#that's throwing errors.

######## Energy-range related variables - may need to be changed between runs ##########

band_check = '2-10'					#This determines what energy band the 
							#fluxes and luminosities are calculated in.
							#Options of (all in keV):
							#'0.5-2'
							#'0.5-7'
							#'0.5-8'
							#'2-10'

if band_check == '0.5-2':				#The band and specific energy to use when
	srcflux_band = '0.5:2:1.56' #soft		#calculating fluxes. It's of the form
elif band_check == '0.5-7':				#(lower):(upper):(specific_energy), with
	srcflux_band = '0.5:7:2.3' #broad		#everything in keV. Since the specific
elif band_check == '0.5-8':				#energy changes depending on which band is
	srcflux_band = '0.5:8:2.3' #~broad		#used, it seemed easier to just do if 
elif band_check == '2-10':				#statements then letting it be user-driven.
	srcflux_band = '2:10:4' #hard			#If you add a new energy band, you'll need
else:							#to update this part as well.
	print("band_check variable not recognized")
	sys.exit()


filter_check = '0.5-7'					#Filters the adjusted and cleaned x-ray
							#image that we use to find the final sources
							#and get flux and luminosity and all that.
							#Options of (all in keV):
							#'0.5-2'
							#'0.5-7'
							#'2-7'


######## Galaxy-specific variables - will need to be changed with each new galaxy #############
#r50_all_gals 		= [7.5, 15, 15, 15, 15]		#Petrosian 50% light radius for each galaxy,
							#in arcsec. Ordered from lowest obsid to 
							#highest obsid.

#galdist 		= [17.53, 20, 20, 20, 20]	#Distance to each galaxy (if known) in Mpc,
							#ordered from lowest obsid to highest obsid.

#galdist_flag 		= True				#Whether or not to use the galaxy distances
							#to compute luminosities. True to compute,
							#False to not.




######## Other variables/parameters - these probably won't need much (if any) changing ###########
coord_frame 		= 'icrs'			#What coordinate frame to use when making the
							#galaxy region. I believe Chandra uses 'icrs',
							#which is the default for this. 'fk5' would 
							#likely make a fine choice as well, the
							#differences would be very, very small.

gal_exclude_rad 	= 3				#Multiplicative factor - this * r50 (the 
							#Petrosian 50% light radius) defines the 
							#galaxy region, used to separate out sources
							#of interest.

fluxim_astrom_psfecf 	= 0.9				#For the astrometry fluximage, make a psf with
							#an enclosed energy fraction of this much.
							#Using 0.9 helps to avoid finding faint
							#sources.

wavdet_astrom_scales 	= '1.0 1.4 2.0 2.8 4.0'		#Scales to use for the astrometry wavdetect.
							#Pretty sure these are in pixels. For more 
							#info, can check the wavdetect ciao thread.

catalog_filter 		= '$imag<22'			#How to filter the catalog sources retrieved
							#from SDSS. The format is like
							#'$column_name<value' (or >). Things like
							#imag or rmag (for column name) generally work
							#though you can check the actual column names
							#via the catalog tool in ds9.

astrom_match_radius 	= 3				#For astrometry, in arcseconds - we require 
							#the separation between a wavdetect and a 
							#catalog source to be less than this in order
							#to count as a match. (can also play around
							#with this a bit)

astrom_match_residlim 	= 2				#For astrometry, in arcseconds - this is the
							#maximum residual we can have (distance
							#between already matched sources) in order to
							#keep the sources as a good match.

astrom_match_method 	= 'trans'			#How to get the astrometrical corrections - 
							#options are 'rst' (for a full rotation,
							#scale factor, and translation) or 'trans'
							#(for just a translation). Since the
							#astrometry corrections seem to tend to be
							#very small (sub-pixel), keeping this as
							#'trans' works best.

min_matches 		= 1				#Minimum number of matches (between wavdetect
							#and catalog sources) required in order to 
							#correct the astrometry. Should be at least
							#1 (obviously).

bin_time_var 		= 200				#How to bin along the time axis when finding
							#background flares. 200 is what is used in
							#the ciao thread, and seems to work well.

fluximage_psfecf 	= 0.393				#For analysis - make a psf with this enclosed
							#energy fraction.

wavdet_analysis_scales 	= '1.0 1.4 2.0 2.8 4.0'		#Scales to use for analysis wavdetect. Units
							#are, again, pixels.

analysis_psfecf 	= 0.9				#For making the source regions (used to
							#measure flux), which are regions
							#corresponding to a psf with this enclosed
							#energy fraction, (cont. below)

analysis_energy 	= 4.5				#at this specific energy (in keV).

bg_radius_factor 	= 12				#Multiplicative factor - background regions
							#(for the source regions above) are annuli
							#with inner radius of the source region and
							#outer radius of this factor * the radius of
							#the source region. Should make pretty big
							#(10-12?) to get good coverage of the galaxy.

srcflux_PhoIndex	 = '1.8'			#For flux - we use a power law spectral model
							#with a photon index of this variable.

minflux_counts 		= 2				#Input 'source counts' used to calculate the
							#expected minimum detectable flux for each
							#galaxy.

minflux_bgcounts 	= 0				#Input 'background counts' used to calculate
							#the expected minimum detectable flux for
							#each galaxy.

ds9_reg_color 		= 'magenta'			#What color to make the output ds9 source and
							#corresponding positional error regions.

round_dec 		= 2 				#How many decimals to round to when printing
							#out data (applies to flux values in a
							#strange way - read on)

flux_ref 		= 1e15				#When printing out flux values, we print out
							#(original flux value)*(this value) and round
							#that to round_dec number of decimal places.

mult_hdrs 		= False				#All sources (from all the galaxies) are 
							#combined into one file (with appropriate
							#galaxy and obsid noted) - this is whether
							#to reprint the column headers (counts, flux,
							#etc.) each time a new galaxy's worth of 
							#sources is added to the file.


#Checking that the band_check and filter_check aren't horribly mismatched
if filter_check == '0.5-2':
	if band_check != '0.5-2':
		print('Filter and band mismatch - please fix')
		print('Filter: '+filter_check)
		print('Band:   '+band_check)
		sys.exit()
if filter_check == '2-7':
	if band_check != '2-10':
		print('Filter and band mismatch - please fix')
		print('Filter: '+filter_check)
		print('Band:   '+band_check)
		sys.exit()

#We initialize a number of lists which correspond to galaxy-specific information, like whether 
#the astrometry was corrected, and if so by how much, and how many matches, etc.
astrom_corrected_list = 	list([]) 
all_gal_name_list = 		list([])
all_gal_ra_list = 		list([])
all_gal_dec_list = 		list([])
a11_list = 			list([])
a12_list = 			list([])
a21_list = 			list([])
a22_list = 			list([])
t1_list = 			list([])
t2_list = 			list([])
gal_astrom_matches = 		list([])
gal_any_xray_srcs = 		list([])
gal_min_flux_list = 		list([])
gal_exp_bg_srcs_list = 		list([])
gal_good_srcs_list = 		list([])
gal_bad_srcs_list = 		list([])
src_summary_df_list = 		list([])
src_all_df_list = 		list([])


#We also get the names of the various directories, for use later. We'll need to return to the code
#directory to find the error on the counts, so we get that. Next, we go up out of the code directory
#into the directory with all the obsids. We save a list of the obsid directories, and save that
#parent directory to return to later. Additionally, we make the directory we'll save the results in
#now, for convenience.
code_directory = os.getcwd()
results_folder = '/results_'+band_check+'-band__'+filter_check+'-filter'
sp.run(["mkdir", "-p", code_directory + results_folder + "/regions_and_images"]) #for later use
sp.run(["mkdir", "-p", code_directory + results_folder + "/cleaned_lightcurves"]) #for later use

os.chdir('..') #going up a directory
all_obsids = next(os.walk('.'))[1] #get list of obsids
all_obsids.sort()
all_obsids = all_obsids[:-1] #as os.walk gets the xray_flux_code dir too, so we delete that entry
original_directory = os.getcwd()


#############################################################
########################### START ###########################
#############################################################

for n in range(0,len(all_obsids)):
	
	os.chdir(original_directory) #bring us back to the original directory
	obsid = all_obsids[n] #focus on one specific galaxy
	os.chdir(obsid) #bring us into that galaxy's obsid directory
	
	#############################################################
	####################### INFO and REPRO ######################
	#############################################################
	'''
	We first find the evt2 file and unzip it if it's zipped. We then get basic info on the galaxy 
	- which is assumed to be the target of the observation - such as the name, RA, and Dec, and
	print it to the screen. Then we reprocess the data with chandra_repro, making sure to check 
	whether things are in faint or vfaint mode. If we detect that chandra_repro has already been 
	run, then we skip it - no need to waste time doing it again. After that we move into the
	repro directory to get started on astrometry matching.
	'''
	
	#Finding the evt2 file. If there's an unzipped evt2 file, we use it. If there isn't, we check
	#to see if there's a zipped file, unzip it, and use it. If neither of these two are found we
	#throw an error and halt. The glob command returns a list of all filenames matching the search
	#parameters.
	test_evt2_unzipped = glob.glob('primary/*evt2.fits')
	test_evt2_zipped = glob.glob('primary/*evt2.fits.gz')
	
	if len(test_evt2_unzipped) == 1:
		evt2_file = test_evt2_unzipped[0]
	elif len(test_evt2_zipped) == 1:
		gunzip(test_evt2_zipped[0])
		evt2_file = test_evt2_zipped[0][:-3] #to get rid of the '.gz'
	else:
		sys.exit('evt2 file not found!')
	
	
	#Now we'll get target galaxy's coordinates and name, which we'll get from the header. Note 
	#that this assumes that the galaxy is the target of the observation. 
	data,hdr = fits.getdata(evt2_file, header=True)
	galname = hdr['OBJECT']
	galra = hdr['RA_TARG']
	galdec= hdr['DEC_TARG']
	galcoords = SkyCoord(ra=galra*u.degree, dec=galdec*u.degree, frame = coord_frame)
	all_gal_name_list.append(galname)
	all_gal_ra_list.append(galra)
	all_gal_dec_list.append(galdec)
	
	
	#Printing the info to the screen. We print the RA and dec in both hms and degrees.
	#We also print the filter_check and band_check information.
	print('--------------------------------------------------------')
	print(obsid,' - ',galname)
	print(galcoords.to_string('hmsdms',precision=1))
	print(galra,galdec)
	print('Filter on final X-ray image:  '+filter_check+' keV')
	print('Energy band for finding flux: '+band_check+' keV')
	
	
	#Now we check whether we need to reprocess the data, and do it if we need to. We also make
	#sure to find whether to use FAINT or VFAINT mode, and store this information for later.
	repro_folder = 'repro_'+band_check+'-band__'+filter_check+'-filter'
	chandra_repro.punlearn()
	chandra_repro.indir		= os.getcwd()
	chandra_repro.outdir		= os.getcwd() + '/' + repro_folder
	chandra_repro.badpixel		= True
	chandra_repro.process_events	= True
	chandra_repro.destreak		= True
	chandra_repro.set_ardlib	= True
	
	datamode = dmkeypar(evt2_file,'datamode',echo=True)
	if datamode == 'VFAINT':
		chandra_repro.check_vf_pha = True
		datamode_str = 'VFAINT'
	else:
		chandra_repro.check_vf_pha = False
		datamode_str = 'FAINT'
		
	chandra_repro.pix_adj		= 'default' #default is EDSER for ACIS (see thread)
	chandra_repro.cleanup		= False     #keep intermediary files
	chandra_repro.clobber		= False
	
	if repro_folder not in os.listdir():
		print('Datamode is ' + datamode)
		print('Reprocessing data...',end=' ',flush=True)
		chandra_repro()
		print('Done')
	else:
		print('Data already reprocessed - skipping chandra_repro')
	
	
	#And now move into the repro directory.
	os.chdir(repro_folder)
	
	#Define bpixfile and set it in ardlib. Probably not necessary, as adding this doesn't seem
	#to actually change any of the results.
	bpixfile = 'acisf' + obsid +'_repro_bpix1.fits'
	ardlib.punlearn()
	acis_set_ardlib(bpixfile)
	
	
	
	#############################################################
	#################### ASTROMETRY MATCHING ####################
	#############################################################
	'''
	**Note - we focus on only the s3 chip from here on out**
	We want to match sources in our image to sources from the SDSS catalog to improve the
	astrometry. We'll do this by using wavdetect on our image to find sources then get the 
	catalog sources from SDSS, then match the two, excluding any wavdetect sources within the  
	galaxy region. If we get enough matches (see VARIABLES section to change the minimum number 
	of matches required -  min_matches) then we update the astrometry of our evt2 file. 
	Otherwise, we copy the old evt2file over without updating the astrometry. We take note of how 
	many matches we found and how much we adjusted the astrometry by, for later.
	'''
	
	#We first get the filenames of the evt2 file and the fov file - the latter is for focusing
	#on just the s3 chip. We'll later end up copying these files and changing the new files to
	#exclude the galaxy - we'll go ahead and get/store those filename strings now, for
	#convenience.
	evt2_repro = glob.glob('acis*_repro_evt2.fits')[0]
	
	evt2_repro_fov = glob.glob('acis*_repro_fov1.fits')[0]
	evt2_repro_fov_s3 = evt2_repro_fov[:-5] + '_s3.fits'
	
	#We now need to make the galaxy region that we'll use to exclude wavdetect sources.
	#This involves mostly messing around with formatting in order to get the ciao functions to
	#work. We first get the ra and dec in celestial coordinates (hms), then turn it into colon
	#format (hh:mm:ss). We then get the radius of circular region in arcminutes. Finally, we
	#make the region string in a format that ciao will work with - e.g. for a circular region,
	#circle(02:30:42.5, 05:20:31.5, 0.5') where the ' denotes arcmin - and make the region. We 
	#also get the area of the region (in pixels) for later use.
	cel_ra = galcoords.ra.to_string(u.hour)
	cel_dec = galcoords.dec.to_string(u.degree,alwayssign=True)
	
	cel_ra = cel_ra.replace('h',':').replace('m',':').replace('s','')
	cel_dec = cel_dec.replace('d',':').replace('m',':').replace('s','')
	
	r50 = r50_all_gals[n] #in arcseconds
	radius = gal_exclude_rad*r50/60 #to get in arcminutes.
	
	ciao_region = 'circle(' + cel_ra + ',' + cel_dec +',' + str(radius)+'\')'
	
	dmmakereg.punlearn()
	dmmakereg.region 	= ciao_region
	dmmakereg.outfile	= 'gal_exclude_phys.reg'
	dmmakereg.kernel	= 'ascii'
	dmmakereg.wcsfile	= evt2_repro
	dmmakereg.clobber	= True
	dmmakereg()
	
	gal_reg = region.CXCRegion('gal_exclude_phys.reg') #using the region package
	gal_area = gal_reg.area()
	
	
	#Note that while *technically* this region should work fine in DS9, sometimes it can be a bit
	#off, not quite sure why. So we also make another region file which is the exact same region,
	#just in a format that DS9 will play nice with.
	wcs_region = reg.CircleSkyRegion(center=galcoords,radius = radius*u.arcmin)
	wcs_region.write('gal_ds9_region.reg', format='ds9', overwrite=True)
	
	
	#Here we make successive copies of the evt2_repro file so that the end file includes only
	#the s3 chip. Could probably be done in one command, but whatever.
	dmcopy.punlearn()
	dmcopy.infile = evt2_repro_fov + '[ccd_id=7]'
	dmcopy.outfile = evt2_repro_fov_s3
	dmcopy.clobber = True
	dmcopy() #Including only s3 chip - step 1
	
	dmcopy.punlearn()
	dmcopy.infile = evt2_repro + '[sky=region(' + evt2_repro_fov_s3 + ')]'
	dmcopy.outfile = evt2_repro[:-5] + '_s3.fits'
	dmcopy.clobber = True
	dmcopy() #Including only s3 chip - step 2

	
	#Now we run first fluximage and then wavdetect. We keep the fluximage bands set (i.e. not a
	#variable) and use the psfecf as specified in the VARIABLES section - fluxim_astrom_psfecf.
	#We check to see if fluximage has been run for the purpose of astrometry-fixing, and run it
	#if it hasn't. We run wavdetect on the result of fluximage, using the scales defined in the
	#VARIABLES section (wavdet_astrom_scales). We check to see if wavdetect has been run for
	#the purposes of fixing the astrometry, and run it if it hasn't.
	fluximage.punlearn()
	fluximage.infile	= evt2_repro[:-5] + '_s3.fits'
	fluximage.outroot	= 'fixastrom_nogal/'
	fluximage.bands		= '0.5:7:2.3' #basically just the broad band
	fluximage.binsize	= 1
	fluximage.psfecf	= fluxim_astrom_psfecf
	fluximage.clobber	= True
	if 'fixastrom_nogal' not in os.listdir():
		print('(Astrometry) Starting fluximage...',end=' ',flush=True)
		fluximage()
		print('Done')
	else:
		print('(Astrometry) Fluximage already run - skipping')
	
	wavdetect.punlearn()
	wavdetect.infile 	= 'fixastrom_nogal/0.5-7_thresh.img'
	wavdetect.psffile	= 'fixastrom_nogal/0.5-7_thresh.psfmap'
	wavdetect.expfile	= 'fixastrom_nogal/0.5-7_thresh.expmap'
	wavdetect.outfile	= 'fixastrom_nogal/wav_astrom_srcs.fits'
	wavdetect.scellfile	= 'fixastrom_nogal/scell.fits'
	wavdetect.imagefile	= 'fixastrom_nogal/imgfile.fits'
	wavdetect.defnbkgfile	= 'fixastrom_nogal/nbgd.fits'
	wavdetect.scales	= wavdet_astrom_scales
	wavdetect.sigthresh	= 1e-06 #running on s3 chip only
	wavdetect.clobber	= True
	wavdet_outfile_check = glob.glob('fixastrom_nogal/wav_astrom_srcs.fits')
	if len(wavdet_outfile_check) == 1:
		print('(Astrometry) Wavdetect already run - skipping')
	else:
		print('(Astrometry) Starting wavdetect...',end=' ',flush=True)
		wavdetect()
		print('Done')
	
	
	#Now excluding the wavdetect sources that are within our galaxy region - we do this source
	#by source, by turning the regions into python objects via the region package to examine
	#their properties. We get the ids of the wavdetect regions inside the galaxy region, and copy
	#the wavdetect source region file over using dmcopy to exclude those regions (the [1:-1] is
	#used to get rid of the brackets from converting the list to a string). We have to wrap that
	#copying part in an if statement because if there are no inside sources (astrom_id_list is
	#empty) then excluding #row= just deletes every source, which is not great.
	astrom_wav_regions = region.CXCRegion('fixastrom_nogal/wav_astrom_srcs.fits')
	gal_astrom_reg     = region.CXCRegion('gal_exclude_phys.reg')
	
	astrom_id_list = list([])
	for m in range(0,len(astrom_wav_regions)):
		#Using .shapes gives us a list of objects, and we pick out one to focus on.
		indiv_reg = astrom_wav_regions.shapes[m]
		
		xp = indiv_reg.xpoints
		yp = indiv_reg.ypoints #getting x,y coords of the center of the wavdetect source
		
		#If the wavdetect region's center is inside the galaxy region, add its id to the list
		if gal_astrom_reg.is_inside(xp,yp):
			astrom_id_list.append(indiv_reg.component)
	
	if len(astrom_id_list) > 0:
		dmcopy.punlearn()
		dmcopy.infile	= 'fixastrom_nogal/wav_astrom_srcs.fits[exclude #row='+str(astrom_id_list)[1:-1]+']'
		dmcopy.outfile	= 'fixastrom_nogal/wav_astrom_out_srcs.fits'
		dmcopy.clobber	= True
		dmcopy()
	else:
		dmcopy.punlearn()
		dmcopy.infile	= 'fixastrom_nogal/wav_astrom_srcs.fits'
		dmcopy.outfile	= 'fixastrom_nogal/wav_astrom_out_srcs.fits'
		dmcopy.clobber	= True
		dmcopy()
	
	#We take the wavdetect sources and filter them to require significance greater than 5 - it's
	#recommended in the thread about fixing absolute astrometry
	# http://cxc.harvard.edu/ciao4.11/threads/reproject_aspect/
	dmcopy.punlearn()
	dmcopy.infile		= 'fixastrom_nogal/wav_astrom_out_srcs.fits[SRC_SIGNIFICANCE=5:]'
	dmcopy.outfile		= 'fixastrom_nogal/wav_astrom_out_srcs_sigGT5.fits'
	dmcopy.clobber		= True
	dmcopy()
	
	
	#Now we get the SDSS catalog sources. We do this via command line interface with DS9, using
	#the same file that we ran wavdetect on. We also check whether we already retrieved the
	#catalog sources or not, and only retrieve them if we need to.
	filename = os.getcwd() + '/fixastrom_nogal/0.5-7_thresh.img'
	cat_file_check = glob.glob('fixastrom_nogal/sdss_cat_srcs.tsv')
	if len(cat_file_check) == 1:
		print('(Astrometry) SDSS catalog sources already retrieved - skipping')
	else:
		print('(Astrometry) Retrieving SDSS catalog sources...',end=' ',flush=True)
		catalog_filter_terminal = "'"+catalog_filter+"'"
		sp.run(["ds9", filename, "-catalog", "sdss", "-catalog", "filter", catalog_filter, "-catalog", "export", "tsv", "fixastrom_nogal/sdss_cat_srcs.tsv", "-exit"])
		print('Done')
	
	
	#And now we find the number of matching sources. While it would be nice to be able to use
	#reproject_aspect for this, that only works if we have 3 or more matching sources, and 
	#sometimes we might only get 1 or 2 good matching sources and we might still like to fix the
	#astrometry based on those. Also, we have to call wcs_match from the command line. This is
	#because if we find no matches (or sometimes if we just find less than three), wcs_match
	#throws an error and, if we use the python-integrated version, will just exit the program.
	#In order to get around this, we open a shell and run it from the command line, suppressing
	#all console output. Then we can just check the logfile, and only correct the wcs if we
	#found more matches than min_matches.
	print('(Astrometry) Matching wavdetect and catalog sources...',end=' ',flush=True)
	wcs_infile	= 'fixastrom_nogal/wav_astrom_out_srcs_sigGT5.fits'
	wcs_refsrcfile	= 'fixastrom_nogal/sdss_cat_srcs.tsv[opt skip=1][cols ra=col3,dec=col4]'
	wcs_outfile	= 'fixastrom_nogal/out.xform'
	wcs_wcsfile	= 'fixastrom_nogal/0.5-7_thresh.img' #keep same file as used for wcs_update
	wcs_logfile	= 'fixastrom_nogal/wcsmatch_log.txt'
	wcs_radius	= str(astrom_match_radius) #in arcsec
	wcs_residlim	= str(astrom_match_residlim) #in arcsec
	wcs_method	= str(astrom_match_method)
	wcs_clobber	= 'yes'
	wcs_verbose	= '5'
	
	
	#The CIAO function wcs_match is throwing out errors (and sometimes segfaults) when
	#there are no matches from the source file to the catalog file (so no matching 
	#sources). Fixed this real quick with the skip_obsid_astrom, where you put in a list
	#of obsids to skip the astrometry on. This works not terribly because the wcs_match
	#errors only happen when there are no matches anyways, so setting the nmatches to
	#zero shouldn't hurt much.
	if obsid not in skip_obsid_astrom:
		wcs_match.punlearn() #just in case
		sp.run(["punlearn", "wcs_match"]) #also just in case
		sp.run(["wcs_match", "infile="+wcs_infile, "refsrcfile="+wcs_refsrcfile, "outfile="+wcs_outfile, "wcsfile="+wcs_wcsfile, "logfile="+wcs_logfile, "radius="+wcs_radius, "residlim="+wcs_residlim, "method="+wcs_method, "clobber="+wcs_clobber, "verbose="+wcs_verbose], stdout=sp.DEVNULL, stderr=sp.DEVNULL) #$***$
		#The use of DEVNULL at the end there suppresses the output, warnings and errors.
		
		#Finding the number of matches and storing these for later use.
		wcs_match_log = open('fixastrom_nogal/wcsmatch_log.txt','r')
		filetext = wcs_match_log.read()
		wcs_match_log.close()
		
		
		#Had trouble since for one galaxy had '-1 matches remaining' and previously had the regex
		#expression without the . in it, which wasn't picking up that minus sign. This way should
		#fix that problem while still working for the normal positive numbers of matches.
		match = re.search('(.\d+) sources remain',filetext)
		nmatches = int(match.group(1))
	else:
		print("Skipping astrometry - galaxy excluded via skip_obsid_astrom")
		nmatches = 0
	print(str(nmatches) + ' matches found')
	gal_astrom_matches.append(nmatches) #saving for later use

	
	#And now we update the astrometry, which means fixing the astrometry in the evt2 fits file,
	#updating the asol files, and updating the asol keyword in the evt2 fits file to point to the
	#new asol files. So first we copy the evt2_repro file over, focusing on the s3 chip.
	evt2_repro_cor = evt2_repro[:-5] + '_corrected.fits'
	
	dmcopy.punlearn()
	dmcopy.infile 		= evt2_repro + '[ccd_id=7]' #again, focusing only on the s3 chip
	dmcopy.outfile 		= evt2_repro_cor
	dmcopy.opt 		= 'all' #to copy all the fits references in the header
	dmcopy.clobber 		= True
	dmcopy()
	
	
	#To get the asol file, we go into the header of the evt2_repro_cor file (which is the
	#file we will be updating the astrometry for) and pull out the name of the file.
	#Originally, we searched for asol files in the repro file directory, but recently
	#there were issues with multiple asol files, only one of which was actually being
	#used in the fits file we care about, so we switched to doing it this way. Note that
	#this will work only if there is a single asol file being used for the observation.
	#It'd be fairly easy to update to work with multiple asol files, but I don't know what
	#the keywords for that would be so I haven't put that in the code.
	#This way also avoids the issue with updating already updated asol files that the last
	#method ran into. For example; say you'd previously run the code and ended up with an
	#evt2_repro_cor file that had an updated asol file (with the _new string appended). 
	#When running the code again, the code will replace the previous evt2_repro_cor file
	#(that had the corrected asol file as the ASOLFILE keyword) with a new evt2_repro_cor
	#file with the old, uncorrected asol file as the keyword. Then that corrected asol
	#file from the previous run will get replaced in the new run by a new corrected asol
	#file. So there's no duplicates.
	asol_file = dmkeypar(evt2_repro_cor, 'ASOLFILE', echo=True)
	
	
	#Now, we only update the astrometry if we get a number of matches greater than our min_matches
	#variable. Otherwise, we'll just leave everything as is.
	if nmatches >= min_matches:
		print('Correcting astrometry...',end=' ',flush=True)
		#We correct the asol file - by making a new one with updated astrometry.
		wcs_update.punlearn()
		wcs_update.infile	= asol_file
		wcs_update.outfile	= asol_file[:-5] + '_new.fits'
		wcs_update.transformfile= 'fixastrom_nogal/out.xform'
		wcs_update.wcsfile	= 'fixastrom_nogal/0.5-7_thresh.img'#same as wcs_match
		wcs_update.logfile	= 'fixastrom_nogal/wcsupdate_asol'+str(m+1)+'_log.txt'
		wcs_update.clobber	= True
		wcs_update.verbose	= 5
		wcs_update()
		
		#And now correcting the fits file. Well, really what we do is update the copy of the
		#file we previously made.
		wcs_update.punlearn()
		wcs_update.infile 	= evt2_repro_cor
		wcs_update.outfile	= '' #as we just update the file, not create a new file
		wcs_update.transformfile= 'fixastrom_nogal/out.xform'
		wcs_update.wcsfile	= 'fixastrom_nogal/0.5-7_thresh.img'
		wcs_update.logfile	= 'fixastrom_nogal/wcsupdate_evt2_log.txt'
		wcs_update.clobber	= True
		wcs_update.verbose	= 5
		wcs_update()
		
		#Now updating the ASOLFILE keyword in the .fits file to point to the corrected asol
		#file.
		dmhedit.punlearn()
		dmhedit.infile		= evt2_repro_cor
		dmhedit.filelist	= 'none'
		dmhedit.operation	= 'add'
		dmhedit.key		= 'ASOLFILE'
		dmhedit.value		= asol_file[:-5] + '_new.fits'
		dmhedit()
		
		astrom_corrected_list.append(1) #making note that we corrected the astrometry
		
		#We also find how much the astrometry was adjusted by, storing the information for
		#later use. The adjustment is of the form 
		#[xnew,ynew] = [ [a11,a12],[a21,a22] ]*[xold,yold] + [t1,t2]
		#i.e. just like a normal transformation matrix. The dmlist command we use returns one
		#long string that looks like eg. '#a11     a12\n    1.0    0.0' for just a11 and a12,
		#where the method used (astrom_match_method) is 'trans' for just a translation. What
		#we do is split that into a list with two entries - the column labels, and the data. 
		#Then we take the entry that's the data and split it up into an array, turning the 
		#number strings into floats. Then we assign everything out.
		#NOTE all these are in units of pixels.
		dmlist.punlearn()
		transform = dmlist('fixastrom_nogal/out.xform[cols a11,a12,a21,a22,t1,t2]', opt='data,clean')
		transform = str.split(transform,'\n')
		transform = str.split(transform[1])
		transform = np.asarray(transform,float)
		a11 = transform[0]
		a12 = transform[1]
		a21 = transform[2]
		a22 = transform[3]
		t1  = transform[4]
		t2  = transform[5]
		
		a11_list.append(a11)
		a12_list.append(a12)
		a21_list.append(a21)
		a22_list.append(a22)
		t1_list.append(t1)
		t2_list.append(t2)
		
		print('Done')
	
	else:
		#Since we didn't correct the astrometry, we append a 0 to astrom_corrected_list and
		#append zeros to all the astrometry adjustment lists.
		astrom_corrected_list.append(0)
		a11_list.append(0)
		a12_list.append(0)
		a21_list.append(0)
		a22_list.append(0)
		t1_list.append(0)
		t2_list.append(0)
		
		
		print('Not enough matches found for astrometric correction - proceeding without.')
		
	
	
		
	
	
	#############################################################
	################### FILTERING and CLEANING ##################
	#############################################################
	'''
	Here we take our corrected evt2 image and filter it for background flares. This involves 
	excluding the galaxy region and any sources we found e.g. the wavdetect sources. While those
	sources were found before the astrometry correction, the astrometry corrections are generally
	very small (sub-pixel) so they should be fine to use here, as all we're using them for is
	to screen out sources. We also do the standard grade filter (on the full image, before
	excluding any regions) and bin the image, before running deflare on it.
	'''
	
	
	#We first need an updated galaxy region (now that we've corrected the astrometry). We can
	#still use the ciao_region string, as that was pulled from the targeting header keywords,
	#which are independent of the astrometry of the image.
	dmmakereg.punlearn()
	dmmakereg.region	= ciao_region
	dmmakereg.outfile	= 'gal_corrected_phys.fits'
	dmmakereg.kernel	= 'ascii'
	dmmakereg.wcsfile	= evt2_repro_cor #updated to the astrometrically corrected file
	dmmakereg.clobber	= True
	dmmakereg()
	
	
	#Next up is the filtering. The [EVENTS] filter gives us just the events (photons hitting the
	#plate). The grade filter is the standard grade filter - see the linked file in the 'grade'
	#entry in the ciao dictionary on the website, right after Table 6.6. Requiring status=0 gives
	#us only the good events - see the dictionary entry for status. The cols filter as written
	#removes the phas column, which may not be strictly necessary, but whatever, it works.
	dmcopy.punlearn()
	dmcopy.infile	= evt2_repro_cor + '[EVENTS][grade=0,2,3,4,6,status=0][cols -phas]'
	dmcopy.outfile	= 'filtered_full_evt2.fits'
	dmcopy.clobber	= True
	dmcopy()
	
	
	#Next we exclude the galaxy and wavdetect regions, and bin the data.
	dmcopy.punlearn()
	dmcopy.infile	= 'filtered_full_evt2.fits[exclude sky=region(gal_corrected_phys.fits)]'
	dmcopy.outfile	= 'lc_nogal_evt2.fits'
	dmcopy.clobber	= True
	dmcopy()
	
	wavdet_astrom_srcs = 'fixastrom_nogal/wav_astrom_srcs.fits'
	
	dmextract.punlearn()
	dmextract.infile	= 'lc_nogal_evt2.fits[exclude sky=region('+wavdet_astrom_srcs+')][bin time=::'+str(bin_time_var)+']'
	dmextract.outfile	= 'lc_binned.fits'
	dmextract.opt		= 'ltc1' #the option to use if you're binning on time, which we are
	dmextract.clobber	= True
	dmextract()
	
	
	#Now we run deflare, and make the final cleaned evt2 file.
	deflare.punlearn()
	deflare.infile	= 'lc_binned.fits'
	deflare.outfile	= 'lc_clean.gti'
	deflare.method	= 'sigma'
	deflare.save	= 'lc_clean_'+ obsid
	deflare()
	
	dmcopy.punlearn()
	dmcopy.infile	= 'filtered_full_evt2.fits[@lc_clean.gti]'
	dmcopy.outfile	= 'acis_full_clean_evt2.fits'
	dmcopy.clobber	= True
	dmcopy()
	
	
	#############################################################
	######################### ANALYSIS ##########################
	#############################################################
	'''
	We filter the image as according to the filter_check user-input variable.
	We then run fluximage and wavdetect, if they haven't already been run for analysis. For 
	fluximage we have a variable controlling the enclosed energy fraction (fluximage_psfecf) 
	and the energy it's at is determined by the filter_check variable. 
	
	We then focus on the wavdetect sources that fall within our galaxy region (note that we count 
	something as within the galaxy region only if its *center* is within the region). We also get 
	the counts in these wavdetect regions for later use in determining the positional error of 
	the sources (see CALCULATING ERRORS section).
	
	We then take these wavdetect sources and make circular regions enclosing a specific energy 
	fraction at a specific energy (e.g. ecf of 90% at 4.5 keV). These are variables (see 
	analysis_psfecf and analysis_energy in the VARIABLES section). We'll use these for finding
	fluxes of the sources. While we're at it, we also find what the radii of circular regions
	enclosing 95% of the energy would be, for later use in determining the upper limit for the
	positional error of the sources (see CALCULATING ERRORS section).
	
	We also find the minimum detectable flux that we would expect to see in the galaxy, and the
	expected number of background sources to fall within our galaxy region.
	
	We then take these source regions, find background regions (annuli around each source, which
	also exclude any other source regions - see bg_radius_factor in VARIABLES section). We plug
	these into srcflux, which has a variety of parameters associated with it - see srcflux_bands
	and srcflux_PhoIndex. We then run srcflux (if it hasn't already been run) and store the
	results. (note we use the Dickey&Lockman absorption column maps)
	'''

	
	#Filtering the image as according to filter_check
	dmcopy.punlearn()
	
	if filter_check == '0.5-2':
		dmcopy.infile	= 'acis_full_clean_evt2.fits[energy=500:2000]'
	elif filter_check == '0.5-7':
		dmcopy.infile	= 'acis_full_clean_evt2.fits[energy=500:7000]'
	elif filter_check == '2-7':	
		dmcopy.infile 	= 'acis_full_clean_evt2.fits[energy=2000:7000]'
	else:
		print("filter_check variable not recognized")
		sys.exit()
	
	dmcopy.outfile	= 'new_evt2.fits'
	dmcopy.clobber	= True
	dmcopy()
	
	
	#Running fluximage
	fluximage.punlearn()
	fluximage.infile	= 'new_evt2.fits'
	fluximage.outroot	= 'analysis/'
	
	if filter_check == '0.5-2':
		fluximage.bands = '0.5:2:1.56'
	if filter_check == '0.5-7':
		fluximage.bands	= '0.5:7:2.3'
	elif filter_check == '2-7':
		fluximage.bands	= '2:7:4'
	else:
		print("filter_check variable not recognized")
		sys.exit()
	
	fluximage.binsize	= 1
	fluximage.psfecf	= fluximage_psfecf
	fluximage.clobber	= True
	if 'analysis' not in os.listdir():
		print('(Analysis) Starting fluximage...',end=' ',flush=True)
		fluximage()
		print('Done')
	else:
		print('(Analysis) Fluximage already run - skipping')
	
	
	#Running wavdetect
	wavdetect.punlearn()
	wavdetect.infile	= 'analysis/'+filter_check+'_thresh.img'
	wavdetect.psffile	= 'analysis/'+filter_check+'_thresh.psfmap'
	wavdetect.expfile	= 'analysis/'+filter_check+'_thresh.expmap'
	wavdetect.outfile	= 'analysis/wav_srcs.fits'
	wavdetect.scellfile	= 'analysis/scell.fits'
	wavdetect.imagefile	= 'analysis/imgfile.fits'
	wavdetect.defnbkgfile	= 'analysis/nbgd.fits'
	wavdetect.scales	= wavdet_analysis_scales
	wavdetect.sigthresh	= 1e-06 #as we're just running on the s3 chip
	wavdetect.clobber	= True
	wavdet_outfile_check_analysis = glob.glob('analysis/wav_srcs.fits')
	if len(wavdet_outfile_check_analysis) == 1:
		print('(Analysis) Wavdetect already run - skipping')
	else:
		print('(Analysis) Starting wavdetect...',end=' ',flush=True)
		wavdetect()
		print('Done')
	
	
	
	#Now finding out which of the wavdetect sources is within our galaxy region. This involves
	#turning the wavdetect regions and the (corrected) galaxy region into python objects via the
	#region package to examine their properties. We get the ids of the wavdetect regions inside
	#the galaxy region, and copy the wavdetect source region file over using dmcopy to keep
	#only those regions (the [1:-1] is used to get rid of the brackets from converting the list
	#to a string).
	wav_regions = region.CXCRegion('analysis/wav_srcs.fits')
	gal_cor_reg = region.CXCRegion('gal_corrected_phys.fits')
	
	id_list = list([])
	for m in range(0,len(wav_regions)):
		#Using .shapes gives us a list of objects, and we pick out one to focus on.
		indiv_reg = wav_regions.shapes[m]
		
		xp = indiv_reg.xpoints
		yp = indiv_reg.ypoints #getting x,y coords of the center of the wavdetect source
		
		#If the wavdetect region's center is inside the galaxy region, add its id to the list
		if gal_cor_reg.is_inside(xp,yp):
			id_list.append(indiv_reg.component)
	
	dmcopy.punlearn()
	dmcopy.infile	= 'analysis/wav_srcs.fits[#row='+str(id_list)[1:-1]+']' 
	dmcopy.outfile	= 'analysis/wav_srcs_inside.fits'
	dmcopy.clobber	= True
	dmcopy()
	
	
	
	#------------------ small detour to get galaxy-specific values like minimum flux
	#Note that some stuff here isn't explained in detail very well - this is because it's 
	#explained later on in the code, when we find the flux for the actual xray sources.
	
	#We do this here because if we detect no X-ray sources within the galaxy, then we can just
	#skip any more analysis and go on to the next galaxy. However, we still want to find some
	#galaxy-specific values before moving on, which is why this block of code is here instead
	#of later in the script.
	
	#And now finding the expected number of background sources to fall within the galaxy region.
	#We'll do this by using Equation 2 from Moretti et al 2003, which gives us the expected
	#number of background sources per square degree. We'll then just multiply this by the area
	#of our galaxy region in degrees. The equation basically takes in a flux, and tells you how
	#many expected background sources with that flux or higher you would expect to see. So what
	#we'll do is calculate our minimum flux - the smallest thing we would measure and count as a
	#source.
	#We would like to be able to find the minimum flux for a galaxy even if it doesn't happen
	#to have any interesting sources in it (e.g. no detected X-ray sources). So what we do
	#is make a source region (at the same energy and ecf as what we would use for flux) at the
	#center of the galaxy. We'll also define a background region (see further on in the code
	#for more explanation as to how this is done - after all, we have to do it again for the
	#actual sources we've found).
	#With that done, we actually just run srcflux on this source. This gives us all sorts of
	#useful information, such as areas of source and background regions, exposure times, and most
	#importantly, a conversion factor that will take us from a count rate to a flux. So we store
	#all of this data, and we run aprates (plugging in the minflux counts and bg_counts) to get
	#a net rate (i.e. what the net rate would be if we had a source located at the center of the
	#galaxy with counts=minflux_counts and bg_counts=minflux_bg_counts). We then extract the
	#net rate and multiply it by the flux conversion factor to get the flux from this hypothetical
	#source.
	psfsize_srcs.punlearn()
	psfsize_srcs.infile	= 'new_evt2.fits' #used for wcs related stuff
	psfsize_srcs.pos	= cel_ra + ' ' + cel_dec #galaxy position in celestial coords
	psfsize_srcs.outfile	= 'analysis/minflux_init_reg.fits'
	psfsize_srcs.energy	= analysis_energy
	psfsize_srcs.ecf	= analysis_psfecf
	psfsize_srcs.clobber	= True
	psfsize_srcs()
	
	roi.punlearn()
	roi.infile	= 'analysis/minflux_init_reg.fits'
	roi.outsrc	= 'analysis/minflux_src.fits'
	roi.radiusmode	= 'mul'
	roi.bkgradius	= bg_radius_factor
	roi.group	= 'exclude'
	roi.targetbkg	= 'target'
	roi.clobber	= 'True'
	roi()
	
	#And making minflux_fin_reg.src.reg and minflux_fin_reg.bg.reg files for use in srcflux
	sp.run(["splitroi", "analysis/minflux_src.fits", "analysis/minflux_fin_reg"])
	
	#Running srcflux - see further on in the code for explanation of stuff
	srcflux.punlearn()
	srcflux.infile		= 'new_evt2.fits' #already filtered
	srcflux.pos		= 'analysis/minflux_init_reg.fits' #for positions only
	srcflux.outroot		= 'minflux_calc/'
	srcflux.bands		= srcflux_band
	srcflux.srcreg		= '@-analysis/minflux_fin_reg.src.reg'
	srcflux.bkgreg		= '@-analysis/minflux_fin_reg.bg.reg'
	srcflux.psfmethod	= 'arfcorr'
	srcflux.conf		= 0.9 #confidence level of 90% (for flux uncertainties)
	srcflux.model		= 'xspowerlaw.pow1'
	srcflux.paramvals	= 'pow1.PhoIndex=' + str(srcflux_PhoIndex)
	srcflux.absmodel	= 'xsphabs.abs1'
	srcflux.absparams	= 'abs1.nh=%GAL%' #using %GAL makes it retriev D&L maps
	srcflux.clobber		= True
	if 'minflux_calc' not in os.listdir():
		print('(Analysis) Starting minimum flux srcflux...',end=' ',flush=True)
		srcflux()
		print('Done')
	else:
		print('(Analysis) Minimum flux srcflux already run - skipping')
	
	#Extracting the things we need for aprates from the flux file and assigning them to variables
	#Note that sometimes using dmlist here would throw warnings out (ignorable, thankfully) which
	#messed up the line splitting thing I did here initially - now we specifically pull out the 
	#last line of the dmlist result, which should be the list of data we want
	dmlist.punlearn()
	minflux_dmlist = dmlist('minflux_calc/_'+band_check+'.flux[cols area,bg_area,psffrac,bg_psffrac, exposure,bg_exposure,umflux_cnv]', opt='data,clean')
	minflux_dmlist = minflux_dmlist.split('\n')
	#minflux_dmlist = minflux_dmlist[1]
	minflux_dmlist = minflux_dmlist[-1] #changed from 1 to -1 to fix warnings problem
	minflux_dmlist = minflux_dmlist.split()
	minflux_area		= float(minflux_dmlist[0])
	minflux_bg_area		= float(minflux_dmlist[1])
	minflux_psffrac		= float(minflux_dmlist[2])
	minflux_bg_psffrac	= float(minflux_dmlist[3])
	minflux_exposure	= float(minflux_dmlist[4])
	minflux_bg_exposure	= float(minflux_dmlist[5])
	minflux_umflux_cnv	= float(minflux_dmlist[6])
	
	
	#Plugging stuff into aprates
	aprates.punlearn()
	aprates.n	= minflux_counts #user variable
	aprates.m	= minflux_bgcounts #user variable
	aprates.A_s	= minflux_area
	aprates.A_b	= minflux_bg_area
	aprates.alpha	= minflux_psffrac
	aprates.beta	= minflux_bg_psffrac
	aprates.T_s	= minflux_exposure
	aprates.T_b	= minflux_bg_exposure
	aprates.outfile	= 'minflux_calc/aprates_netrate_minflux_src.par'
	aprates.clobber	= True
	aprates()
	
	#getting source rate
	tmp_file = open('tmp','w')
	sp.run(["pget", "minflux_calc/aprates_netrate_minflux_src.par", "src_rate"], stdout=tmp_file)
	tmp_file.close()
	temp_minflux_netrate = open('tmp','r').read()
	temp_minflux_netrate = float(temp_minflux_netrate)
	
	#getting minimum flux
	minflux_value = minflux_umflux_cnv * temp_minflux_netrate
	gal_min_flux_list.append(minflux_value) #saving for later use
	
	#Now we plug this into Equation 2 from Moretti et al 2003. They did two fits, for the soft
	#and hard xray bands. We use a variable to choose between them, and if we're not in either
	#of those two bands, we just append nan instead.
	if band_check == '0.5-2' or band_check == '2-10':
		if band_check == '0.5-2':
			a_1	= 1.82
			a_2	= 0.60
			S_band	= 1.48e-14
			N_band	= 6150
		else:
			a_1 	= 1.57
			a_2	= 0.44
			S_band	= 4.5e-15
			N_band	= 5300
			
		N_s = N_band*(2e-15)**a_1/(minflux_value**a_1 + S_band**(a_1-a_2)*minflux_value**a_2)
	
		gal_area_deg = gal_cor_reg.area()*(0.492)**2 #convert from pixels to arcsec
		gal_area_deg = gal_area_deg * (1/3600)**2 #convert to degrees
	
		exp_bg_srcs = N_s * gal_area_deg
		gal_exp_bg_srcs_list.append(exp_bg_srcs) #saving for later use
	else:
		gal_exp_bg_srcs_list.append(np.nan) #if we're not in the soft or hard band
		
	#--------------------end detour
	
	
	
	#There is also the problem of having no xray sources within the galaxy region. We just 
	#continue to the next iteration of the loop right here, and don't bother dealing with 
	#finding flux for specific xray sources (as there aren't any). We also print to the screen
	#that we're doing this, append zeros to the good and bad source lists, and copy the galaxy
	#region and new_evt2.fits files over (and the cleaned lightcurve images, for completion).
	if len(id_list) == 0:
		print('No sources detected in galaxy region - continuing to next galaxy.')
		gal_good_srcs_list.append(0)
		gal_bad_srcs_list.append(0)
		
		galreg_filepath = os.getcwd()+'/gal_ds9_region.reg'
		result_reg_filepath = code_directory+results_folder+'/regions_and_images/'
		result_lc_filepath = code_directory+results_folder+'/cleaned_lightcurves/'
		
		sp.run(["cp", os.getcwd()+"/new_evt2.fits", result_reg_filepath])
		sp.run(["cp", galreg_filepath, result_reg_filepath])
		sp.run(["cp", os.getcwd()+"/lc_clean_"+obsid+".pdf", result_lc_filepath])
		
		temp_current_dir = os.getcwd() #defining for convenience
		temp_resreg_dir = code_directory+results_folder+'/regions_and_images'
		temp_lc_dir = code_directory+results_folder+'/cleaned_lightcurves'
		temp_galname = galname.replace(' ','_')
		temp_new_evt2_filename = obsid+'_'+temp_galname+'_new_evt2.fits'
		temp_new_galreg_filename = obsid+'_'+temp_galname+'_galaxy_region_ds9.reg'
		temp_new_lc_filename = obsid+'_'+temp_galname+'_clean_lightcurve.pdf'
		
		os.chdir(temp_resreg_dir) #change directory
		sp.run(["mv", "new_evt2.fits", temp_new_evt2_filename])
		sp.run(["mv", "gal_ds9_region.reg", temp_new_galreg_filename])
		os.chdir(temp_lc_dir) #change to lightcurve directory
		sp.run(["mv", "lc_clean_"+obsid+".pdf", temp_new_lc_filename])
		
		os.chdir(temp_current_dir) #change directories back
		continue
	
	
	#Now getting the counts in the wavdetect regions for later use in finding the position errors.
	#Unfortunately, the wavdetect output only lists net counts, and we need original counts, so we
	#have to take the net counts and add it to the background counts, then round to the nearest
	#integer in order to get the counts.
	wavdetect_counts = list([])
	for m in range(0,len(id_list)):
		#First we focus in on a single source, making a new file with only that source in it.
		dmcopy.punlearn()
		dmcopy.infile	= 'analysis/wav_srcs_inside.fits[#row='+str(m+1)+']'
		dmcopy.outfile	= 'analysis/wav_src'+str(m+1)+'_inside.fits'
		dmcopy.clobber	= True
		dmcopy()
		
		#Now we get the net_counts and the background counts, add them up and round to the 
		#nearest integer to get the counts.
		dmlist.punlearn()
		temp_dmlist = dmlist('analysis/wav_src'+str(m+1)+'_inside.fits[cols net_counts,bkg_counts]', opt='data,clean')
		temp_dmlist = temp_dmlist.split()
		temp_wavdet_counts = round( float(temp_dmlist[3]) + float(temp_dmlist[4]) )
		wavdetect_counts.append(temp_wavdet_counts)
	
	
	#Now that that's done with, we make our circular regions enclosing a specific energy fraction
	#at a specific energy (e.g. ecf of 90% at 4.5 keV) that we'll use for analysis.
	psfsize_srcs.punlearn()
	psfsize_srcs.infile	= 'new_evt2.fits' #used for wcs related stuff
	psfsize_srcs.pos	= 'analysis/wav_srcs_inside.fits'
	psfsize_srcs.outfile	= 'analysis/init_psf_srcs.fits'
	psfsize_srcs.energy	= analysis_energy
	psfsize_srcs.ecf	= analysis_psfecf
	psfsize_srcs.clobber	= True
	psfsize_srcs()
	
	
	#And while we're at it we'll also make apertures for 95% ecf as we'll need these radii to
	#put upper limits on the position uncertainty (see ERRORS section). We need these radii in
	#arcseconds, and dmlist gives them in pixels, so we'll need to convert to arcseconds using
	#the conversion of 0.492 pixels to arcseconds that applies to all Chandra observations.
	#See http://cxc.harvard.edu/proposer/POG/html/chap6.html
	#These apertures are at 1.5 keV (see Hong et al 2005 section 2.3)
	psfsize_srcs.punlearn()
	psfsize_srcs.infile	= 'new_evt2.fits' #used for wcs related stuff
	psfsize_srcs.pos	= 'analysis/wav_srcs_inside.fits'
	psfsize_srcs.outfile	= 'analysis/95ecf_psf_srcs.fits'
	psfsize_srcs.energy	= 1.5 #kev
	psfsize_srcs.ecf	= 0.95
	psfsize_srcs.clobber	= True
	psfsize_srcs()
	
	dmlist.punlearn()
	r_psfecf95 = dmlist('analysis/95ecf_psf_srcs.fits[cols r]',opt='data,clean')
	r_psfecf95 = r_psfecf95.split('\n')
	r_psfecf95 = np.array(r_psfecf95[1:],float)
	r_psfecf95 = 0.492*r_psfecf95
	
	
	#Now back to the analysis. We need to turn our source regions into source and background
	#regions, which we'll do using the roi tool. We'll need to separate the source and background
	#regions later, though, as roi puts them all in the same file.
	roi.punlearn()
	roi.infile	= 'analysis/init_psf_srcs.fits'
	roi.outsrc	= 'analysis/indiv_psf_src%02d.fits' #the % numbers the different source files
	roi.radiusmode	= 'mul' #so background radius is found by source radius * multiplication fac.
	roi.bkgradius	= bg_radius_factor #multiplication factor
	roi.group	= 'exclude' #for overlapping sources, overlapping part excluded from src reg
	roi.targetbkg	= 'target' #for overlapping sources, bg region is centered on overlap
	roi.clobber	= True
	roi()
	
	
	#And now we have to split the file roi spat out into separate files for source and background
	#regions. However, for some reason the splitroi function that we want to use for this isn't
	#being recognized by python as an existing function. So instead we'll have to do it directly
	#from the command line, like what we did when getting the catalog sources.
	#This creates two region files - fin_psf_srcs.bg.reg and fin_psf_srcs.src.reg, which is what
	#we'll use in srcflux. Note that the syntax for these region files is not DS9 syntax, so they
	#might not load correctly into DS9 (the source one probably will, but the bg one probably
	#won't, as it'll likely involve excluding overlaps).
	sp.run(["splitroi", "analysis/indiv_psf_src\*fits", "analysis/fin_psf_srcs"])
	
	
	#We now copy the source and background regions over to the results_and_regions folder 
	#for ease of access after the code has been run. We also copy the new_evt2.fits file over.
	#We then change directories to the results_and_regions folder, and rename the region and fits
	#files to include the obsid and galaxy name. We then edit the region files to put them in a
	#ds9-friendly format (though you can still only use them over ciao images, as they're in
	#pixels, not sky coordinates). We also do this for the galaxy region defined (much) earlier.
	#We also move the cleaned lightcurve images over.
	src_psf_filepath = os.getcwd()+'/analysis/fin_psf_srcs.src.reg' #defining for convenience
	bg_psf_filepath = os.getcwd()+'/analysis/fin_psf_srcs.bg.reg'
	galreg_filepath = os.getcwd()+'/gal_ds9_region.reg'
	result_reg_filepath = code_directory+results_folder+'/regions_and_images/'
	result_lc_filepath = code_directory+results_folder+'/cleaned_lightcurves/'
	
	sp.run(["cp", src_psf_filepath, result_reg_filepath]) #copying
	sp.run(["cp", bg_psf_filepath, result_reg_filepath])
	sp.run(["cp", os.getcwd()+"/new_evt2.fits", result_reg_filepath])
	sp.run(["cp", galreg_filepath, result_reg_filepath])
	sp.run(["cp", os.getcwd()+"/lc_clean_"+obsid+".pdf", result_lc_filepath])
	
	temp_current_dir = os.getcwd() #defining for convenience
	temp_resreg_dir = code_directory+results_folder+'/regions_and_images'
	temp_lc_dir = code_directory+results_folder+'/cleaned_lightcurves'
	temp_galname = galname.replace(' ','_')
	temp_new_src_filename = obsid+'_'+temp_galname+'_fin_psf_srcs.src.reg'
	temp_new_bg_filename = obsid+'_'+temp_galname+'_fin_psf_srcs.bg.reg'
	temp_new_evt2_filename = obsid+'_'+temp_galname+'_new_evt2.fits'
	temp_new_galreg_filename = obsid+'_'+temp_galname+'_galaxy_region_ds9.reg'
	temp_new_lc_filename = obsid+'_'+temp_galname+'_clean_lightcurve.pdf'
	
	os.chdir(temp_resreg_dir) #change directory
	
	sp.run(["mv", "fin_psf_srcs.src.reg", temp_new_src_filename]) #rename
	sp.run(["mv", "fin_psf_srcs.bg.reg", temp_new_bg_filename])
	sp.run(["mv", "new_evt2.fits", temp_new_evt2_filename])
	sp.run(["mv", "gal_ds9_region.reg", temp_new_galreg_filename])
	os.chdir(temp_lc_dir) #change to lightcurve directory
	sp.run(["mv", "lc_clean_"+obsid+".pdf", temp_new_lc_filename])
	os.chdir(temp_resreg_dir) #change back to regions directory


	#editing the region files, first bg file, then src file
	with open(temp_new_bg_filename,'r') as f:
		bg_newlines = []
		for line in f.readlines():
			bg_newlines.append(line.split('&')[0])
	with open(temp_new_bg_filename,'w') as f:
		for line in bg_newlines:
			f.write(line+'\n')
	
	with open(temp_new_src_filename,'r') as f:
		src_newlines = []
		for line in f.readlines():
			src_newlines.append(line.split('&')[0])
	with open(temp_new_src_filename,'w') as f:
		for line in src_newlines:
			f.write(line+'\n')
	
	os.chdir(temp_current_dir) #change directories back
	
	
	#Finally, we get the flux. For the model we use to calculate it, we choose a powerlaw and 
	#also specify an absorption model (note that we could put the absorption into the .model
	#parameter, but then we wouldn't get any unabsorbed fluxes - in order to do that, we need
	#to supply the absorption separately - see the ahelp file for srcflux). For the srcreg and
	#bkgreg syntax, see the ciao thread. We pick 'quick' for the psfmethod because it makes it
	#easier to do the aperture correction afterwards. -okay, changed, we pick 'arfcorr'
	#because 'quick' fails when there are overlapping sources. This also necessitates
	#changing the aperture correction to rely on the psffrac retrieved from the srcflux
	#output file, instead of just using the psffrac that was entered as a variable.
	srcflux.punlearn()
	srcflux.infile		= 'new_evt2.fits' #already filtered from 0.5-7 keV
	srcflux.pos		= 'analysis/init_psf_srcs.fits' #for positions only
	srcflux.outroot		= 'flux/'
	srcflux.bands		= srcflux_band
	srcflux.srcreg		= '@-analysis/fin_psf_srcs.src.reg'
	srcflux.bkgreg		= '@-analysis/fin_psf_srcs.bg.reg'
	#srcflux.psfmethod	= 'quick'
	srcflux.psfmethod	= 'arfcorr'
	srcflux.conf		= 0.9 #confidence level of 90% (for flux uncertainties)
	srcflux.model		= 'xspowerlaw.pow1'
	srcflux.paramvals	= 'pow1.PhoIndex=' + str(srcflux_PhoIndex)
	srcflux.absmodel	= 'xsphabs.abs1'
	srcflux.absparams	= 'abs1.nh=%GAL%' #using %GAL makes it retriev D&L maps
	srcflux.clobber		= True
	if 'flux' not in os.listdir():
		print('(Analysis) Starting srcflux...',end=' ',flush=True)
		srcflux()
		print('Done')
	else:
		print('(Analysis) Srcflux already run - skipping')
	
	
	#So now we want to get all the data we need (fluxes, errors, counts, etc) in an array or list,
	#to make it easier to deal with. So first, we have to extract this data, which we do using
	#dmlist.
	#IMPORTANT NOTE - ADDING ANY EXTRA COLUMNS HERE MIGHT BREAK IT, AND CHANGING THE ORDER MIGHT
	#BREAK IT
	#For instance, the 'r' column actually adds two columns, one with the actual radius, and
	#another that's just all zeros. In order to fix this I've put the 'r' variable at the end of
	#the retrieval list so that that column of zeros ends up getting cut off later on and we don't
	#have to worry about it. So putting any new columns after 'r' will mess things up. Also, if a
	#column gets added that introduces an additional column of zeros (like 'r' did) then that will
	#break things too.
	dmlist.punlearn()
	dmlist.infile	= 'flux/_'+band_check+'.flux[cols rapos,decpos,theta,phi,area,exposure,bg_exposure,counts,err_counts,net_counts,net_err, bg_counts,bg_area,psffrac,bg_psffrac,umflux_cnv,net_mflux_aper,net_mflux_aper_lo,net_mflux_aper_hi, net_umflux_aper,net_umflux_aper_lo, net_umflux_aper_hi,r]'
	dmlist.opt	= 'data,clean'
	dmlist.outfile	= 'flux/data_list_init.txt'
	dmlist()
	
	
	#Before we can read any of this in, we have to modify the text file. Because dmlist doesn't
	#export that data in a nice format, like tsv or whatever. We can still organize it into
	#something that python will recognize, though. The problem here is that the very first
	#character of the text file is #, and if we try to read the file in as is that character will
	#count as its very own column, which will mess everything up. So we delete that character via
	#the command line, by copying the original file to new file, excluding that annoying #.
	token = open("flux/data_list_fin.txt", 'w')
	sp.run(["tail", "-c", "+2", "flux/data_list_init.txt"], stdout = token)
	token.close()
	
	#Now we can read in the data. We do it line by line, using .split() to separate the columns.
	token = open('flux/data_list_fin.txt')
	linestoken = token.readlines()
	resulttoken = list([])
	for x in linestoken:
		resulttoken.append(x.split())
	token.close()
	
	
	#We then switch rows and columns to make it easier to deal with, as right now it's kind of
	#annoying. e.g. right now resulttoken is a list of lists, and the first element 
	#resulttoken[0] is a list of all the parameter names ['RAPOS','DECPOS',etc.], and the second
	#element is the values of all the parameters for the first source. We'd like to instead have 
	#each element of resulttoken be a list of the parameter names, and then the parameter value
	#for each source e.g. ['RAPOS', ra of source 1, ra of source 2, etc.]. So we use this method
	#to transpose this matrix thingy. Note that this method has issues with non-square matrices
	#and will cut out some data - and our list of lists (matrix) is non-square due to that extra
	#column of zeros from the 'r' parameter. However, since we called 'r' last in dmlist, the only
	#data that will get cut off is that extraneous column of zeros, so it works out fine.
	data_list = list(map(list, zip(*resulttoken)))
	
	
	#Now we break this big list into lots of smaller lists of the values of individual parameters,
	#so we can deal with them easier. We do this by finding the index of the parameter in the 
	#first element of the original list of lists, resulttoken - remember that the first element
	#was just a list of all the parameter names. This index, when plugged into data_list, will
	#give us the list of the specific parameter name and its values for each source. For the 
	#arrays of the data, we discard the parameter name, and store it in a separate list that
	#we'll use when writing all the data to a text file.
	param_list = resulttoken[0] #gets the list of ['net_counts','net_err',etc.]
	
	ra_ind 		  = param_list.index('RAPOS')
	dec_ind 	  = param_list.index('DECPOS')
	theta_ind 	  = param_list.index('THETA')
	phi_ind 	  = param_list.index('PHI')
	area_ind          = param_list.index('AREA')
	exposure_ind	  = param_list.index('EXPOSURE')
	bg_exposure_ind	  = param_list.index('BG_EXPOSURE')
	counts_ind	  = param_list.index('COUNTS')
	err_counts_ind	  = param_list.index('ERR_COUNTS')
	net_counts_ind    = param_list.index('NET_COUNTS')
	net_err_ind       = param_list.index('NET_ERR')
	bg_counts_ind  	  = param_list.index('BG_COUNTS')
	bg_area_ind       = param_list.index('BG_AREA')
	psffrac_ind       = param_list.index('PSFFRAC')
	bg_psffrac_ind    = param_list.index('BG_PSFFRAC')
	umflux_cnv_ind 	  = param_list.index('UMFLUX_CNV')
	net_mflux_ind     = param_list.index('NET_MFLUX_APER')
	net_mflux_hi_ind  = param_list.index('NET_MFLUX_APER_HI')
	net_mflux_low_ind = param_list.index('NET_MFLUX_APER_LO')
	net_umflux_ind    = param_list.index('NET_UMFLUX_APER')
	net_umflux_hi_ind = param_list.index('NET_UMFLUX_APER_HI')
	net_umflux_low_ind= param_list.index('NET_UMFLUX_APER_LO')
	r_ind 		  = param_list.index('R[2]')
	
	ra 		  = np.asarray(data_list[ra_ind][1:],float)
	dec 		  = np.asarray(data_list[dec_ind][1:],float)
	theta 		  = np.asarray(data_list[theta_ind][1:],float)
	phi 		  = np.asarray(data_list[phi_ind][1:],float)
	area 	   	  = np.asarray(data_list[area_ind][1:],float)
	exposure	  = np.asarray(data_list[exposure_ind][1:],float)
	bg_exposure	  = np.asarray(data_list[bg_exposure_ind][1:],float)
	counts		  = np.asarray(data_list[counts_ind][1:],float)
	err_counts	  = np.asarray(data_list[err_counts_ind][1:],float)
	net_counts 	  = np.asarray(data_list[net_counts_ind][1:],float)
	net_err 	  = np.asarray(data_list[net_err_ind][1:],float)
	bg_counts 	  = np.asarray(data_list[bg_counts_ind][1:],float)
	bg_area 	  = np.asarray(data_list[bg_area_ind][1:],float)
	psffrac 	  = np.asarray(data_list[psffrac_ind][1:],float)
	bg_psffrac 	  = np.asarray(data_list[bg_psffrac_ind][1:],float)
	umflux_cnv	  = np.asarray(data_list[umflux_cnv_ind][1:],float)
	net_mflux 	  = np.asarray(data_list[net_mflux_ind][1:],float)
	net_mflux_hi 	  = np.asarray(data_list[net_mflux_hi_ind][1:],float)
	net_mflux_low	  = np.asarray(data_list[net_mflux_low_ind][1:],float)
	net_umflux 	  = np.asarray(data_list[net_umflux_ind][1:],float)
	net_umflux_hi 	  = np.asarray(data_list[net_umflux_hi_ind][1:],float)
	net_umflux_low 	  = np.asarray(data_list[net_umflux_low_ind][1:],float)
	r 		  = np.asarray(data_list[r_ind][1:],float)
	
	ra_str 		  = data_list[ra_ind][0]
	dec_str 	  = data_list[dec_ind][0]
	theta_str 	  = data_list[theta_ind][0]
	phi_str 	  = data_list[phi_ind][0]
	area_str    	  = data_list[area_ind][0]
	exposure_str	  = data_list[exposure_ind][0]
	bg_exposure_str	  = data_list[bg_exposure_ind][0]
	counts_str	  = data_list[counts_ind][0]
	err_counts_str	  = data_list[err_counts_ind][0]
	net_counts_str 	  = data_list[net_counts_ind][0]
	net_err_str	  = data_list[net_err_ind][0]
	bg_counts_str 	  = data_list[bg_counts_ind][0]
	bg_area_str	  = data_list[bg_area_ind][0]
	psffrac_str 	  = data_list[psffrac_ind][0]
	bg_psffrac_str 	  = data_list[bg_psffrac_ind][0]
	umflux_cnv_str	  = data_list[umflux_cnv_ind][0]
	net_mflux_str 	  = data_list[net_mflux_ind][0]
	net_mflux_hi_str  = data_list[net_mflux_hi_ind][0]
	net_mflux_low_str = data_list[net_mflux_low_ind][0]
	net_umflux_str 	  = data_list[net_umflux_ind][0]
	net_umflux_hi_str = data_list[net_umflux_hi_ind][0]
	net_umflux_low_str= data_list[net_umflux_low_ind][0]
	r_str 		  = data_list[r_ind][0]
	
	
	#############################################################
	#################### CALCULATING ERRORS #####################
	#############################################################
	'''
	Here we find the errors on the net counts and source position. We do this via user-made
	functions (made by Rich), and there are some subtleties explained in the comments. These 
	error intervals are 90% confidence levels. While we're doing this we check to see if the 
	sources all pass the significance test (see comments for more explanation as to what the test 
	is. 
	
	We also find the expected number of background sources to fall within our galaxy region.
	
	A few notes regarding the data from srcflux:
	counts 		= original source counts (not background subtracted or aperture corrected)
	net_counts 	= background corrected source counts (but not aperture corrected!)
	count_rate 	= counts/exposure
	net_rate	= net_counts/exposure
	'''
	
	
	#So first we'll calculate the errors on net counts. We do this by plugging in the counts
	#(not background subtracted or aperture corrected) into the user-made error functions
	#kraftcl and gehrelscl that automate finding errors using the methods of Kraft et al. 1991
	#and Gehrels (1986). For sources with counts < 10 we use Kraft and take the background into
	#account (for calculating errors). For sources with counts > 10 we use Gehrels, and assume
	#the background is negligible.
	
	#This gives us the upper and lower limits e.g. counts = 5^{+2.42}_{-1.74} (though note I just
	#made those numbers up, I didn't run them through the error functions). Note that these limits
	#are for the *counts*, not the net_counts, which means they still need to be background 
	#subtracted and aperture corrected, so we do that next.
	
	#For example: if we have counts = 7, bg = 0.5, aperture = 0.9. We plug into kraftcl and get
	#limits of (2.97,11.87) for lower and upper (e.g. 7^{+4.87}_{-4.03}). We then background 
	#subtract and aperture correct everything, counts and limits: (value - 0.5)/0.9. From this
	#we get a value for the net counts of 7.22^{+5.41}_{-4.48}.
	
	
	#Then, we also consider the error bars on the upper and lower intervals, and if the ratio of
	#the two is less than sqrt(2), we consider them to be the same to within uncertainty, so we
	#add them in quadrature and use the result for both upper and lower intervals.
	
	#Following the previous example, the intervals are 5.41 and 4.48 for upper and lower, 
	#respectively. Their ratio is 5.41/4.48 = 1.21 < sqrt(2), so we add them in quadrature to get
	#4.97, which we use for both intervals, giving us a final value of 7.22^{+4.97}_{-4.97}. 
	
	#While we're doing this we also check whether any of the source counts are below the 
	#background to within the 95% confidence level (for Kraft). If for instance, the source 
	#counts are below the background (to within error) we'll mark it as a bad source. We include
	#it in analysis for convenience, but we'll keep track for the final stages of writing
	#everything out to a text document.
	
	bg_counts_in_src = (bg_counts/bg_area*area) #number of bg counts expected in the src aperture
	bg_counts_in_src_cor = bg_counts_in_src/psffrac #corrected for aperture size
	
	net_counts_calc	  = list([]) #used to store bg/aperture corrected net counts
	net_counts_errors = list([]) #used to store bg/aperture corrected net count error intervals
	good_src_list	  = list([]) #1 means good source, 0 means bad source
	
	temp_current_dir = os.getcwd()
	os.chdir(code_directory) #need to be in code dir. for error functions to work properly
	
	for m in range(0,np.size(ra)):
		src_counts 	= counts[m]
		src_bg_counts 	= bg_counts_in_src[m]
		src_psffrac	= psffrac[m]
		
		if src_counts < 10:
			lims = np.array(ps.kraftcl(src_counts,src_bg_counts,1)) #find lims
			
			#Check whether source is above background to with 95% confidence level
			checklims = ps.kraftcl(src_counts,src_bg_counts,2)
			if checklims[0] < src_bg_counts:
				good_src_list.append(0)
			else:
				good_src_list.append(1)
		else:
			lims = np.array(ps.gehrelscl(src_counts,1))
			good_src_list.append(1) #as 95%cl test only matters for kraft
		
		#background/aperture correcting
		cor_lims 	= (lims - src_bg_counts)/src_psffrac
		cor_src_counts	= (src_counts - src_bg_counts)/src_psffrac
		
		#getting upper/lower intervals
		cor_lim_up = cor_lims[1] - cor_src_counts
		cor_lim_lo = cor_src_counts - cor_lims[0]
		
		#checking whether we need to add in quadrature or not
		if max(cor_lim_up/cor_lim_lo,cor_lim_lo/cor_lim_up) < math.sqrt(2):
			single_lim = math.sqrt( (cor_lim_up**2 + cor_lim_lo**2)/2 )
			cor_lim_up = single_lim
			cor_lim_lo = single_lim
		
		#Getting final intervals and adding them (and net counts) to the lists
		fin_lims = np.array([cor_lim_lo,cor_lim_up])
		net_counts_calc.append(cor_src_counts)
		net_counts_errors.append(fin_lims)
	
	os.chdir(temp_current_dir) #move back to normal directory after calculating errors
	
	
	#Now finding the error on the source positions. We do this using the formula in Hong et al.
	# 2005 (equation 5), which depends on both the off-axis distance of the source and the 
	#counts (from the wavdetect regions, specifically). In the paper they also put an upper limit
	#on the positional error of no more than the aperture described by 95% ecf at 1.5 keV. We 
	#check whether we had to use the upper limit and store that for writing out later. Note that 
	#theta is in arcminutes, and the resulting positional uncertainty is in arcseconds.
	pos_error	= list([])
	pos_error_up_lim= list([]) #value of 1 means 95%ecf limit was used, 0 means it wasn't
	
	for m in range(0,np.size(theta)):
		src_theta = theta[m]
		src_wav_counts = wavdetect_counts[m] #from back in the ANALYSIS section
		src_r_psfecf95 = r_psfecf95[m] #from back in the ANALYSIS section
		
		if src_wav_counts > 0:
			cf1 = 1/math.log10(1+src_wav_counts) #counts factor - used in the formula
			cf2 = 1/math.log10(2+src_wav_counts)
			cf3 = 1/math.log10(3+src_wav_counts)
			src_pos_err = 0.25 + 0.1*cf1*(1 + cf1) + 0.03*(cf2*src_theta)**2 + 0.0006*(cf3*src_theta)**4
			if src_pos_err > src_r_psfecf95:
				src_pos_err = src_r_psfecf95
				pos_error_up_lim.append(1)
			else:
				pos_error_up_lim.append(0)
		else:
			src_pos_err = src_r_psfecf95
			pos_error_up_lim.append(1)
		
		pos_error.append(src_pos_err)
		
	
	
	#############################################################
	###################### SAVING RESULTS #######################
	#############################################################
	'''
	Here we write out the information we've gathered so far. We write a summary file concerning
	the sources we found, with things like sky position, net counts, fluxes, etc. - the things
	that we're likely to care the most about. We also write an 'everything' file, which contains
	all the information about the sources we have - i.e. all the parameters we pulled from
	srcflux. Note that for these two files we round things appropriately because otherwise it gets
	difficult to read.
	
	We make a file that summarizes all the stuff relating to the galaxies - if they were
	astrometrically corrected, how many matches there were for the astrometric correction, the 
	astrometric corrections themselves, the minimum flux expected, whether any analysis sources
	were deemed 'bad', and probably other stuff I'm forgetting right now. If there was pertinent
	information that varied from galaxy to galaxy (and not from source to source) then it goes in
	this file.
	
	We make a file that summarizes all the parameters used in the code - like the wavdetect
	scales, the astrometry match radius, etc.
	
	We also make a little folder with the relevant regions in it - the galaxy region, the source
	and background regions, and the final positional source regions (with the position error).
	
	We make another folder with the cleaned lightcurve images in it.
	'''
	
	
	#Going back to the code directory and making the results folder in it if the folder doesn't
	#exist yet, and then going into the results folder. Considering we already make the regions
	#folder back at the start of the code, this probably isn't necessary any more, but the code
	#works right now so I don't want to change it.
	os.chdir(code_directory)
	sp.run(["mkdir", "-p", code_directory + results_folder])
	
	os.chdir(code_directory+results_folder)
	
	
	#We now make the summary file. We'll first just use information from all the sources, and then
	#afterwards exclude any bad sources. We'll do this by putting everything in a pandas dataframe
	#and then turn that dataframe into a string and save it to a file.
	bad_src_ind = np.where(good_src_list == 1)[0] #for later use
	num_good_srcs = len(good_src_list) - len(bad_src_ind) #for later use
	flxref = '(' + str(1/flux_ref) + ')'
	
	#Pre-emptive lists for obsid and galaxy
	ident_list_obsid = [obsid]*np.size(ra)
	ident_list_gal = [galname]*np.size(ra)
	
	#Do it this way to keep the columns in the correct order
	od = collections.OrderedDict()
	od['OBSID']			= ident_list_obsid
	od['GALAXY']			= ident_list_gal
	od['SRC']			= ident_list_gal #placeholder updated later in code
	od[ra_str] 			= ra
	od[dec_str] 			= dec
	od['POS_ERROR'] 		= np.round(pos_error,round_dec)
	od['P_ERR_LIM']			= pos_error_up_lim
	od[net_counts_str] 		= np.round(net_counts_calc,round_dec)
	od['NET_COUNTS_ERR']		= [(round(x[0],round_dec),round(x[1],round_dec)) for x in net_counts_errors] #keeps list format while rounding
	od['BG_COUNTS_IN_SRC'] 		= np.round(bg_counts_in_src_cor,round_dec)
	od[net_umflux_str+flxref] 	= np.round(flux_ref*net_umflux,round_dec)
	if galdist_flag == True:
		lum = net_umflux*4*math.pi*(np.asarray(galdist[n])*3.086e24)**2
		loglum = np.log10(lum)
		od['LOG_LUMINOSITY'] 	= np.round(loglum,round_dec)
	
	summary_df = pd.DataFrame(od) #Turn into dataframe
	summary_df = summary_df.drop(bad_src_ind) #exclude bad sources
	
	#put in the source names e.g. X1, X2, etc.
	pd_src_list 	  = list(np.arange(num_good_srcs) + 1) 
	pd_src_list 	  = ['X'+str(k) for k in pd_src_list]
	summary_df['SRC'] = pd_src_list
	
	#Now we either write the dataframe out to the file as a string (if we want multiple headers)
	#or we append the dataframe to a list, if we want everything under one header.
	if mult_hdrs == True:
		if n == 0:
			summary_df_str = summary_df.to_string(justify='left',index=False)
			summary_file   = open('src_summary.txt','w')
			summary_file.write(summary_df_str)
			summary_file.close()
		else:
			summary_df_str = summary_df.to_string(justify='left',index=False, header=mult_hdrs)
			summary_file   = open('src_summary.txt','a') #append to summary.txt
			summary_file.write('\n\n' + summary_df_str) #start new data on new line
			summary_file.close()
	else:
		src_summary_df_list.append(summary_df)
	
	#We now make the file with all the data in it. In order to get the rounding right, we're
	#unfortunately going to have to do it in the same ordered dictionary approach as used for
	#just the summary, as otherwise the fluxes (which are very, very small) tend to get rounded
	#to zero.
	od_all = collections.OrderedDict()
	od_all['OBSID']			 = ident_list_obsid
	od_all['GALAXY']		 = ident_list_gal
	od_all['SRC']			 = ident_list_gal #placeholder, changed later on
	od_all[ra_str]			 = ra
	od_all[dec_str]			 = dec
	od_all['POS_ERROR']		 = np.round(pos_error,round_dec)
	od_all['P_ERR_LIM']		 = pos_error_up_lim
	od_all[theta_str]		 = np.round(theta,round_dec)
	od_all[phi_str]			 = np.round(phi,round_dec)
	od_all[area_str]		 = np.round(area,round_dec)
	od_all[exposure_str]		 = np.round(exposure,round_dec)
	od_all[bg_exposure_str]		 = np.round(bg_exposure,round_dec)
	od_all[counts_str]		 = counts
	od_all[net_counts_str]		 = np.round(net_counts_calc,round_dec)
	od_all['NET_COUNTS_ERR']	 = [(round(x[0],round_dec),round(x[1],round_dec)) for x in net_counts_errors]
	od_all[bg_counts_str]		 = bg_counts
	od_all[bg_area_str]		 = np.round(bg_area,round_dec)
	od_all['BG_COUNTS_IN_SRC']	 = np.round(bg_counts_in_src_cor,round_dec)
	od_all[psffrac_str]		 = np.round(psffrac,round_dec)
	od_all[bg_psffrac_str]		 = np.round(bg_psffrac,round_dec)
	od_all[umflux_cnv_str]		 = umflux_cnv
	od_all[net_umflux_str+flxref]	 = np.round(flux_ref*net_umflux,round_dec)
	od_all[net_umflux_low_str+flxref]= np.round(flux_ref*net_umflux_low,round_dec)
	od_all[net_umflux_hi_str+flxref] = np.round(flux_ref*net_umflux_hi,round_dec)
	if galdist_flag == True:
		od_all['LOG_LUMINOSITY'] = np.round(loglum,round_dec)
	od_all['R (arcsec)']		 = np.round(0.492*r,round_dec) #to go from pix->arcsec
	
	#Turn into dataframe, exclude bad sources, update 'SRC' column with source names. Also print
	#to the screen how many good and bad sources there are.
	all_df 		= pd.DataFrame(od_all)
	all_df 		= all_df.drop(bad_src_ind)
	all_df['SRC']   = pd_src_list
	print(str(len(pd_src_list))+' good sources, '+str(len(bad_src_ind))+' bad sources')
	gal_good_srcs_list.append(len(pd_src_list)) #saving these two for later use
	gal_bad_srcs_list.append(len(bad_src_ind))
	
	
	#Now we either write the dataframe out to the file as a string (if we want multiple headers)
	#or we append the dataframe to a list, if we want everything under one header.
	if mult_hdrs == True:
		if n == 0:
			all_df_str = all_df.to_string(justify='left',index=False)
			all_file = open('src_all.txt','w')
			all_file.write(all_df_str)
			all_file.close()
		else:
			all_df_str = all_df.to_string(justify='left',index=False,header=mult_hdrs)
			all_file = open('src_all.txt','a')
			all_file.write('\n\n' + all_df_str)
			all_file.close()
	else:
		src_all_df_list.append(all_df)
	
	
	#Now making DS9 regions of the sources (point regions) and their positional errors
	#(circular regions). Note that we reference good_src_list, and only make the region if the
	#source counts as good (i.e. above the background level to within the 95% confidence level).
	os.chdir(temp_resreg_dir)
	ds9_srcs = list([])
	k = 0 #used to label the regions - defined separately from m in case there's a bad source
	for m in range(0,np.size(good_src_list)):
		#Check whether the source is good or not - have to do first because of the label.
		check_good = good_src_list[m]
		if check_good == 1:
			k = k+1
		
		#First the point regions
		center_coords = SkyCoord(ra[m],dec[m],unit='deg',frame='fk5')
		ds9_pt_reg = reg.PointSkyRegion(center=center_coords)
		ds9_pt_reg.visual['color'] = ds9_reg_color
		ds9_pt_reg.visual['marker'] = '+'
		
		#Now the positional error regions
		ds9_er_reg = reg.CircleSkyRegion(center=center_coords,radius=pos_error[m]*u.arcsec)
		ds9_er_reg.visual['color'] = ds9_reg_color
		ds9_er_reg.visual['dash']  = '1'
		ds9_er_reg.meta['text']   = 'X' + str(k)
		
		#And store them in the list if the source is good.
		if check_good == 1:
			ds9_srcs.append(ds9_pt_reg)
			ds9_srcs.append(ds9_er_reg)
		
	#ds9_srcs is now a list of region objects, which we need to turn into a region object which
	#is itself a list of regions in order to be able to write it out to a file. Had to update
	#this one to account for the new syntax in the regions package as well. We also make a
	#variable of the name, for cleanliness.
	ds9_srcs_out = reg.Regions(ds9_srcs)
	ds9_srcs_filename = obsid+'_'+temp_galname+'_xray_ds9_regions.reg'
	ds9_srcs_out.write(ds9_srcs_filename, format='ds9', overwrite=True)


#We first cd into the proper directory.
os.chdir(code_directory+results_folder)


#If we didn't want multiple headers then here we stick all the dataframes in the list together into
#one big dataframe, and print that out to the file as a string.
if mult_hdrs == False and len(src_all_df_list) > 0:
	all_df_tot = pd.concat(src_all_df_list)
	all_df_tot_str = all_df_tot.to_string(justify='left',index=False)
	all_file = open('src_all.txt','w')
	all_file.write(all_df_tot_str)
	all_file.close()

	summary_df_tot = pd.concat(src_summary_df_list)
	summary_df_tot_str = summary_df_tot.to_string(justify='left',index=False)
	summary_file = open('src_summary.txt','w')
	summary_file.write(summary_df_tot_str)
	summary_file.close()
elif len(src_all_df_list) <= 0:
	all_file = open('src_all.txt','w')
	all_file.write("No X-ray sources found!")
	all_file.close()

	summary_file = open('src_summary.txt','w')
	summary_file.write("No X-ray sources found!")
	summary_file.close()

#Quickly defining this flxref variable here. Previously, it was only defined in the for loop 
#which led to problems if there were no xray sources detected in any galaxy (since no detected
#sources meant skipping the rest of the for loop, so this variable never got defined).
flxref = '(' + str(1/flux_ref) + ')'

#Writing the galaxy info file. Note that the astrometry matching parameters a11,a12,a21,a22 only
#matter if you're doing a method other than 'trans', hence the exclusion.
od_gal = collections.OrderedDict()
od_gal['OBSID']		= all_obsids
od_gal['GAL_NAME']	= all_gal_name_list
od_gal['GAL_RA']	= all_gal_ra_list
od_gal['GAL_DEC']	= all_gal_dec_list
od_gal['ASTROM_MATCHES']= gal_astrom_matches
od_gal['ASTROM_COR']	= astrom_corrected_list
if astrom_match_method != 'trans':
	od_gal['a11']	= a11_list
	od_gal['a12']	= a12_list
	od_gal['a21']	= a21_list
	od_gal['a22']	= a22_list
od_gal['t1']		= t1_list
od_gal['t2']		= t2_list
od_gal['MIN_FLUX'+flxref]= np.round(flux_ref*np.asarray(gal_min_flux_list),round_dec)
if galdist_flag == True:
	gal_lum = np.asarray(gal_min_flux_list)*4*math.pi*(np.asarray(galdist)*3.086e24)**2
	gal_loglum = np.log10(gal_lum)
	od_gal['MIN_LOG_LUM'] 	= np.round(gal_loglum,round_dec)
od_gal['EXP_BG_SRCS']	= np.round(gal_exp_bg_srcs_list,4)
od_gal['GOOD_SRCS']	= gal_good_srcs_list
od_gal['BAD_SRCS']	= gal_bad_srcs_list

gal_df		= pd.DataFrame(od_gal)
gal_df_str	= gal_df.to_string(justify='left',index=False)
gal_data_file	= open('galaxy_data.txt','w')
gal_data_file.write(gal_df_str)
gal_data_file.write('\n\n\n\n'+':::ASTROM_MATCHES - number of matches found between wavdetect sources and catalog sources when attempting to correct astrometry. \n:::ASTROM_COR - 1 means astrometry was succesfully corrected, 0 means there weren\'t enough matches to correct the astrometry and everything was left as-is. \n:::a11,a12,a21,a22,t1,t2 - the parameters regarding how the astrometry was corrected. The general form is [xnew,ynew] = [ [a11,a12],[a21,a22] ]*[xold,yold] + [t1,t2], where everything is in units of PIXELS. If the astrom_match_method variable was set to \'trans\' (meaning translation), then we only print t1,t2 to the file, as then the rotation matrix consisting of a11,a12,a21,a22 is trivial. \n:::MIN_FLUX - the minimum flux a source could have where we would still see it - this is calculated from the minflux_counts and minflux_bgcounts parameters. \n:::MIN_LOG_LUM - luminosity of minimum flux - only shows up if galaxy distances have been input and galdist_flag has  been set to True. \n:::EXP_BG_SRCS - the expected number of background sources in the galaxy region, as calculated from Moretti et al. 2003. \n:::GOOD_SRCS - the number of good sources in the galaxy - to qualify as \'good\' a source has to be above the background level to within the 95% confidence level, using the methods of Kraft et al. 1991.')
gal_data_file.close()


#We now go through and add comments to the src_all and src_summary text files, to explain what each
#of the columns are for.
summary_file = open('src_summary.txt','a')
summary_file.write('\n\n\n\n'+':::POS_ERROR - positional error on the sources, in arcsec, from Hong et al 2005. (eqn 5). \n:::P_ERR_LIM - Hong et al. 2005 use an upper limit on the positional error of no more than the aperture described by 95% ecf (at 1.5 keV). A 1 here means that the upper limit was used, a 0 means that the error was smaller than the upper limit. \n:::NET_COUNTS - the net counts, background and aperture corrected (note that the NET_COUNTS output from srcflux is only background corrected, not aperture corrected - we use the fully corrected one here, even though we still call it NET_COUNTS). \n:::NET_COUNTS_ERR - this one is actually pretty complicated. The gist of it is we calculate the errors based off of 90% confidence levels using the methods of Kraft et al. 1991 (for counts<10) and Gehrels 1986 (for counts>10). This gives us upper and lower error intervals, e.g. if NET_COUNTS is 3.01 and NET_COUNTS_ERR is (2.61,4.05), this means the net counts value (with error) is 3.01^{+4.05}_{-2.61}. See comments in the code for more detail (CALCULATING ERRORS section). \n:::BG_COUNTS_IN_SRC - number of background counts expected within each source extraction circle, corrected for aperture size. \n:::NET_UMFLUX_APER - the unabsorbed flux calculated using srcflux, in cgs units (erg/s/cm**2) in multiples of the number in parenthesis after it. e.g. if you have NET_UMFLUX_APER(1e-15) as the column name, and a value of 3.62, then that means that the actual flux is 3.62e-15 erg/s/cm**2. \n:::LOG_LUMINOSITY - this only shows up if you have the galdist_flag variable set to True, and requires manually inputting the distances for each galaxy in the code. Just the log luminosity in erg/s.')
summary_file.close()

all_file = open('src_all.txt','a')
all_file.write('\n\n\n\n'+'Most of these are pulled straight out of the srcflux output file - see http://cxc.harvard.edu/ciao4.11/ahelp/srcflux.html for definitions. The ones that aren\'t or were changed are detailed below. \n:::POS_ERROR - positional error on the sources in arcsec - see src_summary.txt for more detail. \n:::P_ERR_LIM - whether the upper limit was used for the positional error - see src_summary.txt for more detail. \n:::NET_COUNTS - note that this is NOT from the NET_COUNTS output of srcflux, as that is only background subtracted, and not aperture corrected. The data displayed here is both background and aperture corrected (see src_summary.txt). \n:::NET_COUNTS_ERR - (lower,upper) error intervals on NET_COUNTS - see src_summary.txt for more detail. \n:::BG_COUNTS_IN_SRC - number of background counts expected within each source extraction circle, corrected for aperture size. \n:::NET_UMFLUX_APER - see src_summary.txt; same also applies to the lower and upper bounds on umflux. \n:::LOG_LUMINOSITY - this only shows up if you have the galdist_flag variable set to True, and requires manually inputting the distances for each galaxy in the code. Just the log luminosity in erg/s. \n:::R (arcsec) - radii of source region (originally in pixels, converted to arcsec here).')
all_file.close()


#And now finally writing out the parameter file, with some brief explanations as to the purpose of
#each parameter, and a little warning at the end concerning what happens if you rerun the code after
#changing something important but not deleting the relevant directory (like changing the wavdetect
#scales for astrometry and not deleting the astrometry directory).
flux_ref_str_temp = '%.' + str(round_dec) + 'E' #getting flux_ref in scientific notation
flux_ref_str = flux_ref_str_temp % Decimal(str(flux_ref))

param_init_file = open('params.txt','w')
param_init_file.write('Parameters used for this run (for more in depth explanations of each parameter, see the code):')
param_init_file.close()
param_file = open('params.txt','a')
param_file.write('\n' + 'obsids: ' + str(all_obsids))
param_file.write('\n' + 'r50_all_gals: ' + str(r50_all_gals)+' [Petrosian 50% light radius in arcsec]')
param_file.write('\n' + 'galdist: '+str(galdist)+' [Distance to galaxies in Mpc]')
param_file.write('\n' + 'galdist_flag: '+str(galdist_flag)+' [whether to use distances to find luminosities]')
param_file.write('\n' + 'band_check: '+band_check+' [band (in keV) used for flux measurements]')
param_file.write('\n' + 'srcflux_band: '+srcflux_band+' [band and specific energy (in keV) used for calculating fluxes via srcflux]')
param_file.write('\n' + 'filter_check: '+filter_check+' [energy range (in keV) to filter the cleaned, astrometry-corrected final x-ray image]')
param_file.write('\n' + 'coord_frame: '+coord_frame)
param_file.write('\n' + 'gal_exclude_rad: '+str(gal_exclude_rad)+'[ this * r50 defines the galaxy region]')

param_file.write('\n\n' + 'Astrometry-related parameters:')
param_file.write('\n' + 'fluxim_astrom_psfec: '+str(fluxim_astrom_psfecf)+' [make psf with this fraction of energy enclosed]')
param_file.write('\n' + 'wavdet_astrom_scales: '+wavdet_astrom_scales)
param_file.write('\n' + 'catalog_filter: '+catalog_filter)
param_file.write('\n' + 'astrom_match_radius: '+str(astrom_match_radius)+' [in arcsec]')
param_file.write('\n' + 'astrom_match_residlim: '+str(astrom_match_residlim)+' [in arcsec]')
param_file.write('\n' + 'astrom_match_method: '+astrom_match_method)
param_file.write('\n' + 'min_matches: '+str(min_matches)+' [minimum number of matches to the catalog required to change the astrometry]')

param_file.write('\n\n' + 'Cleaning/Filtering-related parameters:')
param_file.write('\n' + 'bin_time_var: '+str(bin_time_var))

param_file.write('\n\n' + 'Analysis-related parameters:')
param_file.write('\n' + 'fluximage_psfecf: '+str(fluximage_psfecf)+' [make psf with this enclosed energy fraction]')
param_file.write('\n' + 'wavdet_analysis_scales: '+wavdet_analysis_scales)
param_file.write('\n' + 'analysis_psfecf: '+str(analysis_psfecf)+' [used to make the source apertures for measuring flux from the sources. Make an aperture enclosing this energy fraction,]')
param_file.write('\n' + 'analysis_energy: '+str(analysis_energy)+' [at this specific energy (in keV)]')
param_file.write('\n' + 'bg_radius_factor: '+str(bg_radius_factor)+' [multiply source region radius by this to get background region outer radius]')
param_file.write('\n' + 'srcflux_PhoIndex: '+srcflux_PhoIndex+' [the index to use for the power law (Gamma) for finding the flux]')

param_file.write('\n\n' + 'Error-calculation related parameters:')
param_file.write('\n' + 'minflux_counts: '+str(minflux_counts)+' [source counts used to find minimum detectable flux]')
param_file.write('\n' + 'minflux_bgcounts: '+str(minflux_bgcounts)+' [bg counts used to find minimum detectable flux]')
#param_file.write('\n' + 'minflux_method: '+minflux_method)

param_file.write('\n\n' + 'Printing-results related parameters:')
param_file.write('\n' + 'ds9_reg_color: '+ds9_reg_color+' [color for final source position and positional error regions]')
param_file.write('\n' + 'round_dec: '+str(round_dec)+' [decimal place to round to for displaying data]')
param_file.write('\n' + 'flux_ref: '+flux_ref_str+' [displayed flux is (original flux)*(this number) then rounded to round_dec number of decimal places.]')
param_file.write('\n' + 'mult_hdrs: '+str(mult_hdrs))
param_file.write('\n\n\n\n'+'WARNING: Some of these parameters can be wrong if you run the code once, change the parameter, then run it again without deleting the relevant directory. For example: you choose wavdet_astrom_scales = \'1 2 3 4\' and run the code. Then, after it\'s run and created the astrometry directory, you change wavdet_astrom_scales to \'2 4 6 8\' and rerun the code, WITHOUT deleting the astrometry directory. The code will find the astrometry directory, decide that the astrometry has already been done, and skip doing the astrometrical wavdetect - so *in the code* it\'s not using the new scales. However, in this parameter file, wavdet_astrom_scales will show up as having the new value - even though it wasn\'t used in this particular case. So be careful - if you change anything relating to astrometry, analysis, or srcflux, make sure to delete the appropriate directory to get the code to actually redo things using the new value.')
param_file.close()


#=====================================================================================================
#=====================================================================================================
######################################################################################################
######################################################################################################
######################################################################################################
################################################ END #################################################
######################################################################################################
######################################################################################################
######################################################################################################
