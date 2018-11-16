__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 November 15'
"""
HISTORY
2018-11-15 created by T.Enoto 
"""
import os 
import sys
import yaml 
import time
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

NPHASE_DEFAULT = 250 
NPHASE_LIST_DEFAULT = [25,50,200,250]
YMINYMAX_LIST_DEFAULT = [[0.035,0.08],[1.75e-2,3.1e-2],[4.5e-3,8.0e-3],[3.5e-3,6.5e-3]]
YMIN_NORM_DEFAULT = 3.5e-3
YMAX_NORM_DEFAULT = 6.5e-3
XMIN_ZOOM_DEFAULT = 0.90
XMAX_ZOOM_DEFAULT = 1.03

def get_profile(phase_array,nphase=NPHASE_DEFAULT):
	nevt = len(phase_array)
	profile_count, profile_binedge, profile_patches = plt.hist(phase_array,nphase,range=[0.0,1.0],histtype='step')
	profile_bincenter = 0.5*(profile_binedge[1:]+profile_binedge[:-1])	
	x = np.array(profile_bincenter)
	y = np.array(profile_count)
	xe = np.full(len(x),0.5/float(nphase))
	ye = np.array(np.sqrt(profile_count))
	ny = y / float(nevt)
	nye = ye / float(nevt)
	return x,xe,y,ye,ny,nye

def plot_profile(x,xe,y,ye,outpdf,label='',ylabel='',title='',
		xmin=0.0,xmax=2.0,ymin=None,ymax=None,legend_title=''):
	outdir = os.path.dirname(outpdf)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	plt.clf()
	fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
	plt.errorbar(np.hstack((x,x+1.0)),np.tile(y,2),yerr=np.tile(ye,2),
		marker='',color='k',drawstyle='steps-mid',linewidth=1.0,label=label)
	if ymax != None:
		plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
	plt.title(title)
	legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)
	axes.set_xlim(xmin,xmax)
	if ymin != None and ymax != None:
		axes.set_ylim(ymin,ymax)
	axes.set_xlabel('Pulse Phase')			
	axes.set_ylabel(ylabel)
	axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig(outpdf,dpi=300)

def save_profile_fitsfile(x,xe,y,ye,ny,nye,outfitsfile):
	cols = []
	cols.append(fits.Column(name='PULSE_PHASE',format='1D',array=x))		
	cols.append(fits.Column(name='PHASE_ERROR',format='1D',array=xe))		
	cols.append(fits.Column(name='COUNTS',format='K',array=y,unit='counts'))
	cols.append(fits.Column(name='COUNTS_ERROR',format='1D',array=ye,unit='counts'))	
	cols.append(fits.Column(name='NORM_COUNTS',format='1D',array=ny,unit='au'))
	cols.append(fits.Column(name='NORM_COUNTS_ERROR',format='1D',array=nye,unit='au'))		
	# http://docs.astropy.org/en/stable/io/fits/
	# http://docs.astropy.org/en/stable/io/fits/usage/table.html
	hdu_primary = fits.PrimaryHDU()
	hdu_table = fits.BinTableHDU.from_columns(cols,name='PROFILE')
	hdulist = fits.HDUList([hdu_primary,hdu_table])
	hdulist.writeto(outfitsfile,overwrite=True)	

def plot_two_profiles(x1,x1e,y1,y1e,x2,x2e,y2,y2e,outpdf,
		label1='',label2='',ylabel='',title='',legend_title='',
		xmin=0.0,xmax=2.0,ymin=None,ymax=None):

	outdir = os.path.dirname(outpdf)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	plt.clf()
	fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
	plt.errorbar(np.hstack((x1,x1+1.0)),np.tile(y1,2),xerr=np.tile(x1e,2),yerr=np.tile(y1e,2),
		marker='',color='k',drawstyle='steps-mid',linewidth=1.0,label=label1)
	plt.errorbar(np.hstack((x2,x2+1.0)),np.tile(y2,2),xerr=np.tile(x2e,2),yerr=np.tile(y2e,2),
		marker='',color='r',drawstyle='steps-mid',linewidth=1.0,label=label2)	
	if ymax != None:
		plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
	plt.title(title)
	legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)
	axes.set_xlim(xmin,xmax)
	if ymin != None and ymax != None:
		axes.set_ylim(ymin,ymax)
	axes.set_xlabel('Pulse Phase')			
	axes.set_ylabel(ylabel)
	axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
	plt.subplots_adjust(wspace=0, hspace=0)
	plt.savefig(outpdf,dpi=300)

def generate_profile_fitsfile(all_pulse_phase_array,mpgrp_flag_array,ipgrp_flag_array,
		outfitsfile,nphase_list=[200,250]):
	print('--- method: %s ---' % sys._getframe().f_code.co_name)

	outdir = os.path.dirname(outfitsfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	hdu_list = []
	hdu_list.append(fits.PrimaryHDU())

	cols_list = []
	for nphase in nphase_list:
		cols_list.append([])
		all_phase,all_phase_error,all_count,all_count_error,all_norm,all_norm_error = get_profile(all_pulse_phase_array,nphase)	
		cols_list[-1].append(fits.Column(name='PULSE_PHASE',format='1D',array=all_phase))		
		cols_list[-1].append(fits.Column(name='PHASE_ERROR',format='1D',array=all_phase_error))		

		cols_list[-1].append(fits.Column(name='ALL_COUNTS',format='K',array=all_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='ALL_COUNTS_ERROR',format='1D',array=all_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='ALL_NORM_COUNTS',format='1D',array=all_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='ALL_NORM_COUNTS_ERROR',format='1D',array=all_norm_error,unit='au'))		

		mpgrp_phase,mpgrp_phase_error,mpgrp_count,mpgrp_count_error,mpgrp_norm,mpgrp_norm_error = get_profile(all_pulse_phase_array[mpgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='MPGRP_COUNTS',format='K',array=mpgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='MPGRP_COUNTS_ERROR',format='1D',array=mpgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_COUNTS',format='1D',array=mpgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_COUNTS_ERROR',format='1D',array=mpgrp_norm_error,unit='au'))		

		nonmpgrp_phase,nonmpgrp_phase_error,nonmpgrp_count,nonmpgrp_count_error,nonmpgrp_norm,nonmpgrp_norm_error = get_profile(all_pulse_phase_array[~mpgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='NONMPGRP_COUNTS',format='K',array=nonmpgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='NONMPGRP_COUNTS_ERROR',format='1D',array=nonmpgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='NONMPGRP_NORM_COUNTS',format='1D',array=nonmpgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='NONMPGRP_NORM_COUNTS_ERROR',format='1D',array=nonmpgrp_norm_error,unit='au'))		

		ipgrp_phase,ipgrp_phase_error,ipgrp_count,ipgrp_count_error,ipgrp_norm,ipgrp_norm_error = get_profile(all_pulse_phase_array[ipgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='IPGRP_COUNTS',format='K',array=ipgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='IPGRP_COUNTS_ERROR',format='1D',array=ipgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='IPGRP_NORM_COUNTS',format='1D',array=ipgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='IPGRP_NORM_COUNTS_ERROR',format='1D',array=ipgrp_norm_error,unit='au'))		

		nonipgrp_phase,nonipgrp_phase_error,nonipgrp_count,nonipgrp_count_error,nonipgrp_norm,nonipgrp_norm_error = get_profile(all_pulse_phase_array[~ipgrp_flag_array],nphase)	
		cols_list[-1].append(fits.Column(name='NONIPGRP_COUNTS',format='K',array=nonipgrp_count,unit='counts'))
		cols_list[-1].append(fits.Column(name='NONIPGRP_COUNTS_ERROR',format='1D',array=nonipgrp_count_error,unit='counts'))	
		cols_list[-1].append(fits.Column(name='NONIPGRP_NORM_COUNTS',format='1D',array=nonipgrp_norm,unit='au'))
		cols_list[-1].append(fits.Column(name='NONIPGRP_NORM_COUNTS_ERROR',format='1D',array=nonipgrp_norm_error,unit='au'))		

		mpgrp_norm_sub = mpgrp_norm - nonmpgrp_norm
		mpgrp_norm_sub_error = np.sqrt(mpgrp_norm_error**2 + nonmpgrp_norm_error**2)
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_SUB',format='1D',array=mpgrp_norm_sub,unit='au'))
		cols_list[-1].append(fits.Column(name='MPGRP_NORM_SUB_ERROR',format='1D',array=mpgrp_norm_sub_error,unit='au'))	

		extname = 'PROFILE_N%d' % nphase
		hdu_list.append(fits.BinTableHDU.from_columns(cols_list[-1],name=extname))

		hdu_list[-1].header['NPHASE'] = nphase 
		hdu_list[-1].header['NXALL']   = len(all_pulse_phase_array)
		hdu_list[-1].header['NXMPGRP'] = len(all_pulse_phase_array[mpgrp_flag_array])		
		hdu_list[-1].header['NXIPGRP'] = len(all_pulse_phase_array[ipgrp_flag_array])				
		hdu_list[-1].header['NXNMPGRP'] = len(all_pulse_phase_array[~mpgrp_flag_array])		
		hdu_list[-1].header['NXNIPGRP'] = len(all_pulse_phase_array[~ipgrp_flag_array])		
		hdu_list[-1].header['FR_MPGRP'] = float(len(all_pulse_phase_array[mpgrp_flag_array]))/float(len(all_pulse_phase_array))
		hdu_list[-1].header['FR_IPGRP'] = float(len(all_pulse_phase_array[ipgrp_flag_array]))/float(len(all_pulse_phase_array))
	fits.HDUList(hdu_list).writeto(outfitsfile,overwrite=True)	

def plot_fitsfile(profile_fitsfile,nphase,outpdf,ymin,ymax):
	hdu = fits.open(profile_fitsfile)

	extname = 'PROFILE_N%d' % nphase 
	x1  = hdu[extname].data['PULSE_PHASE']
	x1e = hdu[extname].data['PHASE_ERROR']
	y1  = hdu[extname].data['ALL_NORM_COUNTS']
	y1e = hdu[extname].data['ALL_NORM_COUNTS_ERROR']
	y2  = hdu[extname].data['MPGRP_NORM_COUNTS']
	y2e = hdu[extname].data['MPGRP_NORM_COUNTS_ERROR']

	plot_two_profiles(x1,x1e,y1,y1e,x1,x1e,y2,y2e,outpdf,
		label1='All',label2='MPGRP',ylabel='Normalized counts',title='',legend_title='',
		xmin=0.0,xmax=2.0,ymin=ymin,ymax=ymax)

class XrayEvent():
	def __init__(self,eventfits):
		self.eventfits = eventfits 

		if not os.path.exists(self.eventfits):
			sys.stderr.write('file %s does not exist.\n' % self.eventfits)
			quit()
		print('XrayEvent: {} loaded.'.format(self.eventfits))

		try:
			with fits.open(self.eventfits) as self.hdu:
				self.header  = self.hdu['EVENTS'].header				
				self.data    = self.hdu['EVENTS'].data
				self.columns = self.hdu['EVENTS'].columns
		except OSError as e:
			raise 

		self.set_statistics()

	def set_statistics(self):
		self.numofevt_all = len(self.data['TIME'])
		self.numofevt_mpgrp = len(self.data['MPGRP_FLAG'][self.data['MPGRP_FLAG']])
		self.numofevt_ipgrp = len(self.data['IPGRP_FLAG'][self.data['IPGRP_FLAG']])	
		self.fraction_mpgrp = float(self.numofevt_mpgrp)/float(self.numofevt_all)
		self.fraction_ipgrp = float(self.numofevt_ipgrp)/float(self.numofevt_all)		

	def show_statistics(self):
		print('numofevt_all: {:,}'.format(self.numofevt_all))
		print('numofevt_mpgrp: {:,} ({:.2f}%)'.format(self.numofevt_mpgrp,self.fraction_mpgrp*100.0))		
		print('numofevt_ipgrp: {:,} ({:.2f}%)'.format(self.numofevt_ipgrp,self.fraction_ipgrp*100.0))				

class XrayEventList():
	def __init__(self,xrayevt_list_file,outdir='out/merge'):
		self.xrayevt_list_file = xrayevt_list_file
		self.outdir = outdir
		print('xrayevt_list_file: {}'.format(self.xrayevt_list_file))
		print('outdir: {}'.format(self.outdir))

		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)

	def set_xrayevt_list(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		
		self.xrayevt_array = []
		self.hdu_list = []
		self.time_array_list = []
		self.pulse_phase_array_list = []
		self.mpgrp_flag_array_list = []
		self.ipgrp_flag_array_list = []
		for line in open(self.xrayevt_list_file):
			cols = line.split()
			filepath = cols[0]
			self.xrayevt_array.append(filepath)
			self.hdu_list.append(XrayEvent(filepath))
			self.hdu_list[-1].show_statistics()
			self.time_array_list.append(self.hdu_list[-1].data['TIME'])
			self.pulse_phase_array_list.append(self.hdu_list[-1].data['PULSE_PHASE'])
			self.mpgrp_flag_array_list.append(self.hdu_list[-1].data['MPGRP_FLAG'])			
			self.ipgrp_flag_array_list.append(self.hdu_list[-1].data['IPGRP_FLAG'])						
		print('xrayevt_array: {}'.format(self.xrayevt_array))

	def merge_eventlist(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		
		self.time_array = np.hstack(self.time_array_list)
		self.pulse_phase_array = np.hstack(self.pulse_phase_array_list) 
		self.mpgrp_flag_array = np.hstack(self.mpgrp_flag_array_list)
		self.ipgrp_flag_array = np.hstack(self.ipgrp_flag_array_list)

	def set_statistics(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		
		self.numofevt_all = len(self.time_array)
		self.numofevt_mpgrp = len(self.mpgrp_flag_array[self.mpgrp_flag_array])
		self.numofevt_ipgrp = len(self.ipgrp_flag_array[self.ipgrp_flag_array])	
		self.fraction_mpgrp = float(self.numofevt_mpgrp)/float(self.numofevt_all)
		self.fraction_ipgrp = float(self.numofevt_ipgrp)/float(self.numofevt_all)		

	def show_statistics(self):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)	
		print('numofevt_all: {:,}'.format(self.numofevt_all))
		print('numofevt_mpgrp: {:,} ({:.2f}%)'.format(self.numofevt_mpgrp,self.fraction_mpgrp*100.0))		
		print('numofevt_ipgrp: {:,} ({:.2f}%)'.format(self.numofevt_ipgrp,self.fraction_ipgrp*100.0))				

	def generate_profile(self,pulse_phase_array,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)

		phase,phase_error,count,count_error,norm,norm_error = get_profile(pulse_phase_array,nphase)
		if flag_plot:
			outpdf = '%s/%s_cnt.pdf' % (self.outdir,outbase)
			plot_profile(phase,phase_error,count,count_error,outpdf,label=outbase,
				ylabel='Counts')
			outpdf = '%s/%s_norm.pdf' % (self.outdir,outbase)
			plot_profile(phase,phase_error,norm,norm_error,outpdf,label=outbase,
				ymin=YMIN_NORM_DEFAULT,ymax=YMAX_NORM_DEFAULT,ylabel='Normalized counts')			
		if flag_fitsout:
			outfits = '%s/%s.fits' % (self.outdir,outbase)			
			save_profile_fitsfile(phase,phase_error,count,count_error,norm,norm_error,outfits)			
		return [phase,phase_error,count,count_error,norm,norm_error]

	def generate_profile_all(self,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		self.pls_allsum = self.generate_profile(self.pulse_phase_array,outbase,nphase,flag_plot,flag_fitsout)

	def generate_profile_mpgrp(self,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		self.pls_mpgrp = self.generate_profile(self.pulse_phase_array[self.mpgrp_flag_array],outbase,nphase,flag_plot,flag_fitsout)

	def generate_profile_nonmpgrp(self,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		self.pls_nonmpgrp = self.generate_profile(self.pulse_phase_array[~self.mpgrp_flag_array],outbase,nphase,flag_plot,flag_fitsout)

	def generate_profile_ipgrp(self,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		self.pls_ipgrp = self.generate_profile(self.pulse_phase_array[self.ipgrp_flag_array],outbase,nphase,flag_plot,flag_fitsout)

	def generate_profile_nonipgrp(self,outbase,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)
		self.pls_nonipgrp = self.generate_profile(self.pulse_phase_array[~self.ipgrp_flag_array],outbase,nphase,flag_plot,flag_fitsout)

	def check_excess_mp_relative_to_all(self,outpdf):
		plot_two_profiles(
			self.pls_allsum[0],self.pls_allsum[1],self.pls_allsum[4],self.pls_allsum[5],
			self.pls_mpgrp[0],self.pls_mpgrp[1],self.pls_mpgrp[4],self.pls_mpgrp[5],
			outpdf,ylabel='Normalized counts',
			label1='all',label2='mpgrp',ymin=YMIN_NORM_DEFAULT,ymax=YMAX_NORM_DEFAULT)

		zoom_outpdf = outpdf.replace('.pdf','_zoom.pdf')
		plot_two_profiles(
			self.pls_allsum[0],self.pls_allsum[1],self.pls_allsum[4],self.pls_allsum[5],
			self.pls_mpgrp[0],self.pls_mpgrp[1],self.pls_mpgrp[4],self.pls_mpgrp[5],
			zoom_outpdf,ylabel='Normalized counts',
			label1='all',label2='mpgrp',
			xmin=XMIN_ZOOM_DEFAULT,xmax=XMAX_ZOOM_DEFAULT,
			ymin=YMIN_NORM_DEFAULT,ymax=YMAX_NORM_DEFAULT)

	def generate_profile_fitsfile(self,outfitsfile,nphase_list=NPHASE_LIST_DEFAULT,flag_plot=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)

		generate_profile_fitsfile(self.pulse_phase_array,
			self.mpgrp_flag_array,self.ipgrp_flag_array,
			outfitsfile,nphase_list=nphase_list)

		if flag_plot:
			for i in range(len(nphase_list)):
				outpdf = '%s_n%d.pdf' % (os.path.splitext(outfitsfile)[0],nphase_list[i])
				ymin = YMINYMAX_LIST_DEFAULT[i][0]
				ymax = YMINYMAX_LIST_DEFAULT[i][1]			
				plot_fitsfile(outfitsfile,nphase_list[i],outpdf,ymin=ymin,ymax=ymax)

	def generate_randomized_profile_fitsfile(self,outfitsfile,nphase_list=NPHASE_LIST_DEFAULT,flag_plot=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)		

		np.random.seed(seed=int(time.time()))
		generate_profile_fitsfile(self.pulse_phase_array,
			np.random.permutation(self.mpgrp_flag_array),
			np.random.permutation(self.ipgrp_flag_array),
			outfitsfile,nphase_list=nphase_list)

		if flag_plot:
			for i in range(len(nphase_list)):
				outpdf = '%s_n%d.pdf' % (os.path.splitext(outfitsfile)[0],nphase_list[i])
				ymin = YMINYMAX_LIST_DEFAULT[i][0]
				ymax = YMINYMAX_LIST_DEFAULT[i][1]	
				plot_fitsfile(outfitsfile,nphase_list[i],outpdf,ymin=ymin,ymax=ymax)

	def loop_generate_randomized_profile_fitsfile(self,num_of_loop):
		for i in range(num_of_loop):
			outfitsfile = '%s/rand/rand%03d/random%03d.fits' % (self.outdir,i,i)
			self.generate_randomized_profile_fitsfile(outfitsfile)			

	def run(self):
		self.set_xrayevt_list()
		self.merge_eventlist()
		self.set_statistics()
		self.show_statistics()
		#self.generate_profile_all('def/pls_allsum')
		#self.generate_profile_mpgrp('def/pls_mpgrp')
		#self.generate_profile_nonmpgrp('def/pls_nonmpgrp')		
		#self.generate_profile_ipgrp('def/pls_ipgrp')
		#self.generate_profile_nonipgrp('def/pls_nonipgrp')				
		#self.check_excess_mp_relative_to_all('%s/def/pls_allsum_mpgrp_norm.pdf' % self.outdir)
		outfitsfile = '%s/default/default.fits' % (self.outdir)
		self.generate_profile_fitsfile(outfitsfile)
		self.loop_generate_randomized_profile_fitsfile(500)
		#outfitsfile = '%s/rand/random1.fits' % (self.outdir)
		#self.generate_randomized_profile_fitsfile(outfitsfile)
		#outfitsfile = '%s/rand/random2.fits' % (self.outdir)
		#self.generate_randomized_profile_fitsfile(outfitsfile)

		print('\007 \007 \007')
