#!/usr/bin/env python

import os 
import glob
import pandas as pd 

FILE_CRAB_PULSAR_EPHEMERIS = 'data/radio/original/v181001/crab_pulsar_ephemeris.csv'
PATH_TO_RADIO_DATADIR = 'data/radio/original/v181001'
OUTDIR_ROOT = 'data/radio/fits/v181009'
SKIP_DATAID_LIST = []

df = pd.read_csv(FILE_CRAB_PULSAR_EPHEMERIS,
	delim_whitespace=True,skiprows=[0,1,2],header=0)

for index, row in df.iterrows():
	if row['P0ms'] == '--':
		continue
	if str(row['dataid']) in SKIP_DATAID_LIST:
		continue
	print(row['dataid'])

	input_dir = '%s/%s' % (PATH_TO_RADIO_DATADIR,row['dataid'])
	output_dir = '%s/%s' % (OUTDIR_ROOT,row['dataid'])

	if len(glob.glob('%s/*_orgformat_GTI.txt' % input_dir)) > 0:
		gtifile = glob.glob('%s/*_orgformat_GTI.txt' % input_dir)[0]
	elif len(glob.glob('%s/*_GTI.txt' % input_dir)) > 0:
		gtifile = glob.glob('%s/*_GTI.txt' % input_dir)[0]
	else:
		print("can not find gti file.")
		exit()

	outfits = '%s.fits' % os.path.splitext(os.path.basename(gtifile))[0]
	cmd  = 'giantradiopulse/cli/radio.py convert-gti-txt2fits '
	cmd += '%s ' % gtifile
	cmd += '--outfitsfile %s/%s ' % (output_dir,outfits)
	print(cmd);os.system(cmd)

	for infile in glob.glob('%s/*.txt' % input_dir):
		if not ('GTI' in infile):
			outfits = '%s.fits' % os.path.splitext(os.path.basename(infile))[0]
			cmd  = 'giantradiopulse/cli/radio.py convert-grp-txt2fits '
			cmd += '%s ' % infile
			cmd += '--outfitsfile %s/%s ' % (output_dir,outfits)
			print(cmd);os.system(cmd)	
