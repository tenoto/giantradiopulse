#!/usr/bin/env python

import giantradiopulse.xrayprofile 

#f = giantradiopulse.xrayprofile.XrayProfileFitsfile('ni1013010104_0mpu7_cl_nibary_phase_corr.fits')
#f.write_fitsfile_with_normalized_extensions()
f = giantradiopulse.xrayprofile.XrayProfileFitsfile('ni1013010104_0mpu7_cl_nibary_phase_corr.fits')
#f.plot_compared_pulseprofiles()
f.plot_compared_pulseprofiles(lag=0,xmin=0.90,xmax=1.10,ymin=1.5e-2,ymax=2.7e-2)
exit()

#f.evaluate_peak_significance(f.)

f = giantradiopulse.xrayprofile.XrayProfileFitsfile()
f.generate_fitsfile(
	'out/crabgrp/v181026/2017221/2017221_setup.yaml',
	'out/crabgrp/v181025/2017221/ni1013010104_0mpu7_cl_nibary_phase.evt',
	'out/crabgrp/v181025/2017221/2017221_UsdS_MPGRPlistDE430SNge5.5_updated.fits',
	'out/crabgrp/v181025/2017221/2017221_UsdS_IPGRPlistDE430SNge5.5_updated.fits',
	'out/crabgrp/v181025/2017221/2017221_UsdS_DE430NEW2_orgFormat_GTI.fits',
	lagrange=2,nphase=60)
