#!/usr/bin/env python

import os 
import sys 
import numpy as np
import astropy.io.fits as fits 

import matplotlib as mpl
import matplotlib.pyplot as plt 

mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '14'
mpl.rcParams['mathtext.default'] = 'regular'
mpl.rcParams['xtick.top'] = 'True'
mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.05' #'.05'
mpl.rcParams['axes.ymargin'] = '.05'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

NUM_OF_PHASEBIN = 180
INEVT = 'data/xray/v181009/2017221/ni1013010104_0mpu7_cl_nibary_phase_gti_grp.evt'

def get_profile(pulse_phase_array,nbin_phase):
	num_of_pulse = len(np.unique(pulse_phase_array))
	#num_of_pulse = len(pulse_phase_array)

	repeat_pulse_phase_array = pulse_phase_array + 1.0
	twocycle_pulse_phase_array = np.concatenate([pulse_phase_array,repeat_pulse_phase_array], axis=0)
	hist_count, hist_binedge, hist_patch = plt.hist(twocycle_pulse_phase_array,
		bins=2*nbin_phase,histtype="step")
	hist_bincenter = 0.5 * (hist_binedge[1:] + hist_binedge[:-1])
	x = hist_bincenter
	y = hist_count
	yerr = hist_count**0.5
	ynorm = hist_count / num_of_pulse
	ynorm_err = hist_count**0.5 / num_of_pulse
	return x, y, yerr, ynorm, ynorm_err, num_of_pulse

### 
# START MAIN 
### 
try: 
	with fits.open(INEVT) as hdu:
		pulse_phase = hdu['EVENTS'].data['PULSE_PHASE']
		mod_pulse_number = hdu['EVENTS'].data['MOD_PULSE_NUMBER']
		flag_mpgrp5p5 = hdu['EVENTS'].data['MPGRP5.5']
except OSError as e:
	raise 
print("file %s is loaded." % INEVT)

#print(flag_mpgrp5p5)
#print(~flag_mpgrp5p5)
#print(len(pulse_phase))
#print(len(pulse_phase[flag_mpgrp5p5]))

nongrp_phase, nongrp_count, nongrp_error, nongrp_count_norm, nongrp_error_norm, nongrp_numpulse = get_profile(pulse_phase[~flag_mpgrp5p5],NUM_OF_PHASEBIN)
mpgrp_phase, mpgrp_count, mpgrp_error, mpgrp_count_norm, mpgrp_error_norm, mpgrp_numpulse = get_profile(pulse_phase[flag_mpgrp5p5],NUM_OF_PHASEBIN)
ymin = 0.95 * min(nongrp_count_norm)
ymax = 1.05 * max(nongrp_count_norm)

plt.clf()
fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
plt.errorbar(nongrp_phase,nongrp_count_norm,yerr=nongrp_error_norm,
	marker='',color='k',drawstyle='steps-mid',linewidth=1.0,label='Non-GRP average')
plt.errorbar(mpgrp_phase,mpgrp_count_norm,yerr=mpgrp_error_norm,
	marker='',color='r',drawstyle='steps-mid',linewidth=1.0,label='MP-GRP selected')
plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
axes.legend(loc='upper right',shadow=False)
axes.set_xlim(0.0,2.0)
axes.set_ylim(ymin,ymax)
axes.set_xlabel('Pulse Phase')			
#axes.set_ylabel('Counts / bin / Number of X-ray Events')
axes.set_ylabel('Counts / bin / Number of Pulse Numbers')
axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
#axes2 = axes.twiny()
#date1 = datetime.datetime(1996,2,20)
#date2 = datetime.datetime(2018,9,29)
#axes2.set_xlim(date1,date2)
#axes.vlines(58315,1.8,2.7,colors='k')
#axes.vlines(54530,1.8,2.7,colors='k')
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('profile.pdf',dpi=300)

axes.set_xlim(0.9,1.1)
#axes.set_ylim(77000,150000)
plt.savefig('profile_zoom.pdf',dppi=300)

