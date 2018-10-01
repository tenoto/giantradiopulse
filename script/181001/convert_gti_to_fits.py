#!/usr/bin/env python

import os 
import glob

input_dir = 'data/radio/original/v181001/2017221'
output_dir = 'data/radio/fits/v181001/2017221'

for infile in glob.glob('%s/*.txt' % input_dir):
	if 'GTI' in infile:
		outfits = '%s.fits' % os.path.splitext(os.path.basename(infile))[0]
		cmd  = 'giantradiopulse/cli/radio.py convert_gti_txt2fits '
		cmd += '%s ' % infile
		cmd += '--outfitsfile %s/%s ' % (output_dir,outfits)
		print(cmd);os.system(cmd)
	else:
		outfits = '%s.fits' % os.path.splitext(os.path.basename(infile))[0]
		cmd  = 'giantradiopulse/cli/radio.py convert_grp_txt2fits '
		cmd += '%s ' % infile
		cmd += '--outfitsfile %s/%s ' % (output_dir,outfits)
		print(cmd);os.system(cmd)		