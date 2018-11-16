#!/usr/bin/env python

import glob
import numpy as np 
import astropy.io.fits as fits

import scipy.optimize as opt

import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.rcParams['font.family'] = 'Times New Roman'
mpl.rcParams['font.size'] = '15'
mpl.rcParams['mathtext.default'] = 'regular'
#mpl.rcParams['xtick.top'] = 'True'
#mpl.rcParams['ytick.right'] = 'True'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
#mpl.rcParams['axes.grid'] = 'True'
mpl.rcParams['axes.xmargin'] = '.08' #'.05'
mpl.rcParams['axes.ymargin'] = '.10'
mpl.rcParams['savefig.facecolor'] = 'None'
mpl.rcParams['savefig.edgecolor'] = 'None'
mpl.rcParams['savefig.bbox'] = 'tight'

numof_xrays_mpgrp = []
significance_n50 = []
significance_n25 = []
significance_n250 = []
for i in range(7):
	fitsfile = 'out/crabgrp/v181116_merge/input%02d/default/default.fits' % i
	if len(glob.glob(fitsfile)) == 0:
		continue
	print(fitsfile)
	hdu = fits.open(fitsfile)
	numof_xrays_mpgrp.append(hdu[1].header['NXMPGRP']/1e+6)
	for extnum in range(1,len(hdu)):
		peak_index = np.argmax(hdu[extnum].data['ALL_NORM_COUNTS'])
		peak_enhance = hdu[extnum].data['MPGRP_NORM_SUB'][peak_index]
		peak_enhance_error = hdu[extnum].data['MPGRP_NORM_SUB_ERROR'][peak_index]
		peak_enhance_significance = peak_enhance / peak_enhance_error
		print('%d: %.2e %.2e %.3e' % (extnum,peak_enhance,peak_enhance_error,peak_enhance_significance))
		if hdu[extnum].name == 'PROFILE_N50':
			significance_n50.append(peak_enhance_significance)
		if hdu[extnum].name == 'PROFILE_N25':
			significance_n25.append(peak_enhance_significance)

	peak_index = np.argmax(hdu['PROFILE_N250'].data['ALL_NORM_COUNTS'])
	all_count_peak_sum = sum(hdu['PROFILE_N250'].data['ALL_COUNTS'][peak_index-1:peak_index+2])
	mpgrp_count_peak_sum = sum(hdu['PROFILE_N250'].data['MPGRP_COUNTS'][peak_index-1:peak_index+2])
	all_norm_peak_sum = float(all_count_peak_sum)/ float(hdu['PROFILE_N250'].header['NXALL'])
	mpgrp_norm_peak_sum = float(mpgrp_count_peak_sum) / float(hdu['PROFILE_N250'].header['NXMPGRP'])
	all_norm_peak_sum_error = np.sqrt(float(all_count_peak_sum))/ float(hdu['PROFILE_N250'].header['NXALL'])
	mpgrp_norm_peak_sum_error = np.sqrt(float(mpgrp_count_peak_sum)) / float(hdu['PROFILE_N250'].header['NXMPGRP'])	
	peak_enhance = mpgrp_norm_peak_sum - all_norm_peak_sum 	
	peak_enhance_error = np.sqrt(all_norm_peak_sum_error**2+mpgrp_norm_peak_sum_error**2)
	peak_enhance_significance = peak_enhance / peak_enhance_error
	print(all_norm_peak_sum,mpgrp_norm_peak_sum,peak_enhance,peak_enhance_error,peak_enhance_significance)
	significance_n250.append(peak_enhance_significance)

def func(x, k):
     return k * np.sqrt(x) 

# The actual curve fitting happens here
optimizedParameters, pcov = opt.curve_fit(func, numof_xrays_mpgrp, significance_n250);
# Use the optimized parameters to plot the best fit

plt.clf()
fig, axes = plt.subplots(1,1,figsize=(7.0,6.0))
#plt.plot(numof_xrays_mpgrp,significance_n250,'ro',markersize=8.0,label='0.04 x 3',zorder=4,edgecolor='black')
#plt.plot(numof_xrays_mpgrp,significance_n50,'ys',markersize=8.0,label='0.02',zorder=3)
#plt.plot(numof_xrays_mpgrp,significance_n25,'b^',markersize=8.0,label='0.04',zorder=2)
plt.scatter(numof_xrays_mpgrp,significance_n250,marker='o',s=100.0,facecolor='r',label='0.004 x 3',zorder=4,edgecolor='black')
plt.scatter(numof_xrays_mpgrp,significance_n50,marker='s',s=100.0,facecolor='y',label='0.02',zorder=3,edgecolor='black')
plt.scatter(numof_xrays_mpgrp,significance_n25,marker='^',s=100.0,facecolor='b',label='0.04',zorder=2,edgecolor='black')
x = np.arange(0,10,0.1)
plt.plot(x,func(x, *optimizedParameters),'r',linestyle='--',zorder=1);
axes.set_xlabel(r'Number of X-ray events coincidenced with MP-GRPs (10$^{6}$ photons)')
axes.set_ylabel('Significance of the enhancement')
axes.set_xlim(0.0,10.0)
axes.set_ylim(0.0,6.0)
plt.legend(loc='upper left',title='phase width')
plt.savefig('out/crabgrp/v181116_merge/growth_curve_n50.pdf')




