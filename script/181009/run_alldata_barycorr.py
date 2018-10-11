#!/usr/bin/env python

import os 
import glob 
import pandas as pd 

FILE_CRAB_PULSAR_EPHEMERIS = 'data/radio/original/v181001/crab_pulsar_ephemeris.csv'
PATH_TO_NICER_OBSDIR = '/Users/enoto/work/drbv1/reporitory/heasarc/data/nicer/data/obs'
#SKIP_DATAID_LIST = ['2018003']
SKIP_DATAID_LIST = []
OUTDIR_ROOT = 'data/xray/v181009'
NPHASE = 100

with open(FILE_CRAB_PULSAR_EPHEMERIS) as f:
	for line in f:
		cols = line.split()
		if cols[1] == 'RA':
			RA_J2000 = float(cols[3])
		if cols[1] == 'DEC':
			DEC_J2000 = float(cols[3])
		if cols[1] == 'Planetary':
			Planetary_ephemeris = str(cols[3].replace('DE',''))
print("RA  = {:.9}".format(RA_J2000))
print("DEC = {:.9}".format(DEC_J2000))

df = pd.read_csv(FILE_CRAB_PULSAR_EPHEMERIS,
	delim_whitespace=True,skiprows=[0,1,2],header=0)

for index, row in df.iterrows():
	if row['P0ms'] == '--':
		continue
	if str(row['dataid']) in SKIP_DATAID_LIST:
		continue
	print(row['dataid'])

	outdir = OUTDIR_ROOT+'/{}'.format(row['dataid'])
	cmd = 'rm -rf %s; mkdir -p %s' % (outdir,outdir)
	print(cmd);os.system(cmd)

	obsid = str(row['nicerobsid'])
	path_nicer_obsid_list = glob.glob('%s/*/%s' % (PATH_TO_NICER_OBSDIR,obsid))
	if len(path_nicer_obsid_list) == 0:
		print("no correspodning input data file.")
		continue 
	path_nicer_obsid = path_nicer_obsid_list[0]
	print(path_nicer_obsid)

	cmd  = 'rm -rf tmp_nitimeconv.out;'
	cmd += 'nitimeconv.py %d -f mjd -s tt > tmp_nitimeconv.out' % row['MJD']
	print(cmd);os.system(cmd)
	with open('tmp_nitimeconv.out') as f:
		for line in f:
			cols = line.split()
			if len(cols) == 0:
				continue 
			if cols[0] == 'NICER':
				met_epoch_day = float(cols[7])
	met_epoch = met_epoch_day + float(row['tJPLms'])*1e-3
	cmd  = 'rm -rf tmp_nitimeconv.out;'
	print(cmd);os.system(cmd)

	# This calculation of F2 comes from the C code at the end of the
	# explanatory notes for the Jodrell ephemeris
	# f2 = 2.0*p1*p1/(p0*p0*p0)	
	p0 = float(row['P0ms'])*1e-3
	p1 = float(row['Pdot0'])
	f2 = 2.0*p1*p1/(p0*p0*p0)	

	nibaryevt  = 'ni%s_0mpu7_cl_nibary.evt' % obsid
	niphaseevt = 'ni%s_0mpu7_cl_nibary_phase.evt' % obsid
	niprofile = 'ni%s_0mpu7_cl_nibary_phase_profile.fht' % obsid	
	file_script = outdir+'/run_xraybary_{}.sh'.format(row['dataid'])
	with open(file_script,'w') as f:
		dump = """#!/bin/sh -f
nibarytime.py \
%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz \
%.6f %.6f %s/auxil/ni%s.orb.gz \
--refframe ICRS --ephem JPLEPH.%s \
--outfits %s/%s
""" % (path_nicer_obsid,obsid,
	RA_J2000,DEC_J2000,path_nicer_obsid,obsid,Planetary_ephemeris,
	outdir,nibaryevt)
		f.write(dump)

		dump = """
rm -f %s/%s 
faddphase_nu.py \
%s/%s \
%.7f \
%.13f \
--nudot=%6.5fe-15 \
--nu2dot=%.6e \
--outfits %s/%s
""" % (outdir,niphaseevt,
	outdir,nibaryevt,
	met_epoch, float(row['nu0']),float(row['nudot0_15']),f2,
	outdir,niphaseevt)
		f.write(dump)

		dump = """
rm -f %s/%s 
fplot_pulseprofile.py \
%s/%s \
--outfits %s/%s \
--nbin %d \
--colname PULSE_PHASE \
--title "Crab pulse profile ni%s all (MET)"    
""" % (outdir,niprofile,
	outdir,niphaseevt,
	outdir,niprofile,
	NPHASE,obsid)
		f.write(dump)

	cmd = 'chmod +x %s' % file_script
	print(cmd);os.system(cmd)

	cmd = '%s' % file_script
	print(cmd);os.system(cmd)















