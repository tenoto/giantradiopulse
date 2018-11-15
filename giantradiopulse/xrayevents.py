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
	# http://docs.astropy.org/en/stable/io/fits/
	# http://docs.astropy.org/en/stable/io/fits/usage/table.html
	hdu_primary = fits.PrimaryHDU()
	hdu_table = fits.BinTableHDU.from_columns(cols,name='PROFILE')
	hdulist = fits.HDUList([hdu_primary,hdu_table])
	hdulist.writeto(outfitsfile,overwrite=True)	


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

	def set_xrayevt_list(self):
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
		self.time_array = np.hstack(self.time_array_list)
		self.pulse_phase_array = np.hstack(self.pulse_phase_array_list) 
		sel f.mpgrp_flag_array = np.hstack(self.mpgrp_flag_array_list)
		self.ipgrp_flag_array = np.hstack(self.ipgrp_flag_array_list)

	def set_statistics(self):
		self.numofevt_all = len(self.time_array)
		self.numofevt_mpgrp = len(self.mpgrp_flag_array[self.mpgrp_flag_array])
		self.numofevt_ipgrp = len(self.ipgrp_flag_array[self.ipgrp_flag_array])	
		self.fraction_mpgrp = float(self.numofevt_mpgrp)/float(self.numofevt_all)
		self.fraction_ipgrp = float(self.numofevt_ipgrp)/float(self.numofevt_all)		

	def show_statistics(self):
		print('numofevt_all: {:,}'.format(self.numofevt_all))
		print('numofevt_mpgrp: {:,} ({:.2f}%)'.format(self.numofevt_mpgrp,self.fraction_mpgrp*100.0))		
		print('numofevt_ipgrp: {:,} ({:.2f}%)'.format(self.numofevt_ipgrp,self.fraction_ipgrp*100.0))				

	def generate_profile_all(self,nphase=NPHASE_DEFAULT,flag_plot=True,flag_fitsout=True):
		print('--- method: %s ---' % sys._getframe().f_code.co_name)

		self.plf_all_phase,self.plf_all_phase_error,self.plf_all_count,self.plf_all_count_error,self.plf_all_norm,self.plf_all_norm_error= get_profile(self.pulse_phase_array,nphase)
		if flag_plot:
			plot_profile(self.plf_all_phase,self.plf_all_phase_error,
				self.plf_all_count,self.plf_all_count_error,'test.pdf',label='hoge')
		if flag_fitsout:
			save_profile_fitsfile(
				self.plf_all_phase,self.plf_all_phase_error,
				self.plf_all_count,self.plf_all_count_error,
				self.plf_all_norm,self.plf_all_norm_error,
				'test.fits')

	def run(self):
		self.set_xrayevt_list()
		self.merge_eventlist()
		self.set_statistics()
		self.show_statistics()
		self.generate_profile_all()
