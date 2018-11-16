#!/usr/bin/env python

import glob
import numpy as np 
import astropy.io.fits as fits

from scipy import optimize

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

# Equation for Gaussian
def gauss(x, A, mu, sigma):
	return A * np.exp(-(x - mu)**2.0 / (2 * sigma**2))

def get_siginificance_n250(fitsfile):
	hdu = fits.open(fitsfile)
	extname = 'PROFILE_N250'
	peak_index = np.argmax(hdu[extname].data['ALL_NORM_COUNTS'])
	bottom_index = np.argmin(hdu[extname].data['ALL_NORM_COUNTS'])

	all_count_peak_sum = sum(hdu[extname].data['ALL_COUNTS'][peak_index-1:peak_index+2])
	mpgrp_count_peak_sum = sum(hdu[extname].data['MPGRP_COUNTS'][peak_index-1:peak_index+2])
	all_norm_peak_sum = float(all_count_peak_sum)/ float(hdu[extname].header['NXALL'])
	mpgrp_norm_peak_sum = float(mpgrp_count_peak_sum) / float(hdu[extname].header['NXMPGRP'])
	all_norm_peak_sum_error = np.sqrt(float(all_count_peak_sum))/ float(hdu[extname].header['NXALL'])
	mpgrp_norm_peak_sum_error = np.sqrt(float(mpgrp_count_peak_sum)) / float(hdu[extname].header['NXMPGRP'])	
	peak_enhance = mpgrp_norm_peak_sum - all_norm_peak_sum 	
	peak_enhance_error = np.sqrt(all_norm_peak_sum_error**2+mpgrp_norm_peak_sum_error**2)
	peak_enhance_significance = peak_enhance / peak_enhance_error

	nave = 10
	bottom_array = hdu[extname].data['ALL_COUNTS'][bottom_index-1-nave:bottom_index+nave]
	bottom_value = float(sum(bottom_array))/float(len(bottom_array)) / float(hdu[extname].header['NXALL'])
	peak_enhance_relative = peak_enhance / (all_norm_peak_sum - bottom_value)
	#print(peak_enhance_relative)
	return [all_norm_peak_sum,mpgrp_norm_peak_sum,peak_enhance,peak_enhance_error,peak_enhance_significance,peak_enhance_relative]

default_sig_result = get_siginificance_n250('out/crabgrp/v181116_rand/input07/default/default.fits')
default_enhance_relative = default_sig_result[5]*100.0
default_enhance_error = default_sig_result[3]*100.0
default_enhance_significance = default_sig_result[4]
print(0, default_enhance_relative,default_enhance_error,default_enhance_significance)

x_array = []
y_array = []
ye_array = []
ys_array = []
for fitsfile in glob.glob('out/crabgrp/v181116_rand/input07/rand/*/*fits'):
	sig_result = get_siginificance_n250(fitsfile)
	enhance_relative = sig_result[5]*100.0
	enhance_error = sig_result[3]*100.0
	enhance_significance = sig_result[4]
	print(len(x_array)+1,enhance_relative,enhance_error,enhance_significance)
	x_array.append(len(x_array)+1)
	y_array.append(enhance_relative)
	ye_array.append(enhance_error)
	ys_array.append(enhance_significance)

plt.clf()
fig, axes = plt.subplots(1,1,figsize=(7.0,6.0))
plt.subplot(3, 1, 1)
plt.title('(phase 0.004 x 3 bin, MP)')
plt.errorbar(x_array,y_array,fmt='ko',yerr=ye_array)
plt.errorbar([0],[default_enhance_relative],fmt='ro',yerr=[default_enhance_error])
plt.ylabel('Enhancement (%)')
plt.ylim(-4.0,4.0)
plt.subplot(3, 1, 2)
plt.errorbar(x_array,ys_array,fmt='ko')
plt.errorbar([0],[default_enhance_significance],fmt='ro',yerr=[default_enhance_error])
plt.ylim(-6.0,6.0)
plt.ylabel(r'Significance ($\sigma$)')
plt.xlabel('Number of trials (n=0 real data, n>1 fake)')
#plt.subplot(3, 1, 3)
plt.savefig('out/crabgrp/v181116_rand/crab_rand_show.pdf')

nbin_hist = 50
plt.clf()
fig, axes = plt.subplots(1,1,figsize=(6.0,6.0))
y, x, _ = plt.hist(y_array,nbin_hist,range=[-4.0,4.0],histtype='step')
#param = gauss.fit(y_array)
#print(param)
#x = np.arange(-4.0,4.0,0.01)
#pdf = gauss.pdf(x)
#plt.plot(x,pdf)
x=(x[1:]+x[:-1])/2 # for len(x)==len(y)
popt, pcov = optimize.curve_fit(gauss, x, y)
x_fit = np.linspace(x[0], x[-1], 1000)
y_fit = gauss(x_fit, *popt)
print(popt)
mean = popt[1]
sigma = np.fabs(popt[2])
print(mean,sigma)
label_fit='m=%.3f s=%.3f' % (mean,sigma)
plt.plot(x_fit, y_fit, lw=2, color="k",label=label_fit)
significance = (default_enhance_relative - mean)/sigma
label_default = r'%.3f $\sigma$' % significance
plt.hist([default_enhance_relative],nbin_hist,range=[-4.0,4.0],histtype='stepfilled',color='r',label=label_default)
title = 'n=%d' % len(y_array)
plt.legend(loc='upper left',title=title,fontsize=10.0)
plt.xlabel('Enhancement (%)')
plt.ylabel('Number of siumlated sample')
plt.savefig('out/crabgrp/v181116_rand/crab_rand_hist.pdf')









