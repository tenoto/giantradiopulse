#!/usr/bin/env python

import os 
import sys 
import numpy as np
import astropy.io.fits as fits 

def append_column(xray_evtfile,grp_fitsfile):
	hdu_xray = fits.open(xray_evtfile)
	xray_mod_pulse_number = hdu_xray['EVENTS'].data['MOD_PULSE_NUMBER']
	num_of_xrays = len(xray_mod_pulse_number)
	print(xray_mod_pulse_number.dtype)

	hdu_grp = fits.open(grp_fitsfile)
	grp_mod_pulse_number = hdu_grp['GRP'].data['NSEQpulse']
	num_of_grps = len(grp_mod_pulse_number)
	print(grp_mod_pulse_number.dtype)

	sys.stdout.write('Number of X-rays: %d\n' % num_of_xrays)
	sys.stdout.write('Number of GRPs: %d\n' % num_of_grps)

	xray_flag_isin_grp = np.isin(xray_mod_pulse_number,grp_mod_pulse_number)
	num_of_xrays_in_grp = np.sum(xray_flag_isin_grp==True)
	sys.stdout.write('Number of X-rays within GRPs: %d\n' % num_of_xrays_in_grp)

	new_columns = fits.ColDefs([
		fits.Column(name='MPGRP5.5',format='L',array=xray_flag_isin_grp)
		])

	hdu_xray_events_columns = hdu_xray['EVENTS'].columns

	xray_evtfile_grp = '%s_grp.evt' % os.path.splitext(xray_evtfile)[0]
	cmd = 'rm -f %s' % xray_evtfile_grp
	print(cmd);os.system(cmd)

	hdu_primary = fits.PrimaryHDU()
	hdu_events = fits.BinTableHDU.from_columns(hdu_xray_events_columns+new_columns,name='EVENTS')
	hdulist = fits.HDUList([hdu_primary,hdu_events])
	hdulist.writeto(xray_evtfile_grp)

def run(xray_evtfile,grp_fitsfile,grp_gtifile):
	print("xray_evtfile: %s" % xray_evtfile)
	print("grp_fitsfile: %s" % grp_fitsfile)
	print("grp_gtifile: %s" % grp_gtifile)

	xray_evtfile_gti = '%s_gti.evt' % os.path.splitext(xray_evtfile)[0]
	cmd = 'rm -f %s' % xray_evtfile_gti
	print(cmd);os.system(cmd)

	cmd = 'fselect '
	cmd += '%s+1 ' % xray_evtfile
	cmd += '%s ' % xray_evtfile_gti
	cmd += '"gtifilter(\'%s[1]\',BARY_TIME,\'BARY_START\',\'BARY_STOP\')"' % grp_gtifile
	print(cmd);os.system(cmd)

	append_column(xray_evtfile_gti,grp_fitsfile)

#run(xray_evtfile='data/xray/v181009/2017221/ni1013010104_0mpu7_cl_nibary_phase.evt',
#	grp_fitsfile='data/radio/fits/v181009/2017221/2017221_UsdS_MPGRPlistDE430SNge5.5_updated.fits',
#	grp_gtifile='data/radio/fits/v181009/2017221/2017221_UsdS_DE430NEW2_GTI.fits'
#	)

run(xray_evtfile='data/xray/v181009/2018003/ni1013010123_0mpu7_cl_nibary_phase.evt',
	grp_fitsfile='data/radio/fits/v181009/2018003/2018003_UsdS_MPGRPlistDE430SNge5.5_updated.fits',
	grp_gtifile='data/radio/fits/v181009/2018003/2018003_UsdS_DE430NEW2_orgformat_GTI.fits'
	)



