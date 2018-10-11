#!/usr/bin/env python

import sys 
import numpy as np
import astropy.io.fits as fits 

xray_evtfile = 'data/xray/v181009/2017221/ni1013010104_0mpu7_cl_nibary_phase.evt'
radio_gti = 'data/radio/fits/v181009/2017221/2017221_UsdS_DE430NEW2_GTI.fits'
radio_grpfile = 'data/radio/fits/v181009/2017221/2017221_UsdS_MPGRPlistDE430SNge5.5_updated.fits'

hdu_xray = fits.open(xray_evtfile)
xray_mod_pulse_number = hdu_xray['EVENTS'].data['MOD_PULSE_NUMBER']
num_of_xrays = len(xray_mod_pulse_number)
print(xray_mod_pulse_number.dtype)

hdu_radio = fits.open(radio_grpfile)
grp_mod_pulse_number = hdu_radio['GRP'].data['NSEQpulse']
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

outfitsfile = 'test.evt'
hdu_primary = fits.PrimaryHDU()
hdu_events = fits.BinTableHDU.from_columns(hdu_xray_events_columns+new_columns,name='EVENTS')
hdulist = fits.HDUList([hdu_primary,hdu_events])
hdulist.writeto(outfitsfile)