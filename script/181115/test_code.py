#!/usr/bin/env python

import numpy as np
import astropy.io.fits as fits

import matplotlib.pyplot as plt 
import matplotlib as mpl
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

def get_xray_profile(phase_array,nphase=250):
	profile_count, profile_binedge, profile_patches = plt.hist(
		phase_array,nphase,range=[0.0,1.0],histtype='step')
	profile_bincenter = 0.5*(profile_binedge[1:]+profile_binedge[:-1])
	return np.array(profile_bincenter),np.array(profile_count)

inevt = 'out/crabgrp/v181115/2017314/ni1013010110_0mpu7_cl_nibary_phase_0p3to12p0keV_grpflag.evt'

hdu = fits.open(inevt)
#print(hdu['EVENTS'].columns)
mpgrp_flag = hdu['EVENTS'].data['MPGRP_FLAG'] 
phase_array =  hdu['EVENTS'].data['PULSE_PHASE'] 
#print(len(mpgrp_flag))
#print(len(mpgrp_flag[mpgrp_flag]))
#print(len(mpgrp_flag[~mpgrp_flag]))
#nonmpgrp_profile = get_xray_profile_bincenter_array(phase_array[~mpgrp_flag],nphase=250)
profile_all_phase, profile_all_count = get_xray_profile(phase_array,nphase=250)

print(profile_all_phase)
print(profile_all_count)

xmin = 0.0
xmax = 2.0
outpdf = 'test.pdf'

plt.clf()
fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
plt.errorbar(
	np.hstack((profile_all_phase,profile_all_phase+1.0)),
	np.tile(profile_all_count,2),
	#yerr=np.tile(self.hdu['NONMPGRP_NORM'].data[yerr_colname],2),
	marker='',color='k',drawstyle='steps-mid',linewidth=1.0,
	#label='No-GRP (%.3e +/- %.3e)' % (peak_nongrp,peak_err_nongrp)
	)
#plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
#plt.title(title)
#legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)
axes.set_xlim(xmin,xmax)
#axes.set_ylim(ymin,ymax)
axes.set_xlabel('Pulse Phase')			
axes.set_ylabel('Normalized count rate')
axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig(outpdf,dpi=300)
