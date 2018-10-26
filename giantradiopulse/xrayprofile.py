__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 October 26'
"""
HISTORY
2018-10-26 created by T.Enoto 
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

class XrayProfileFitsfile():
	""" Represents 
	"""

	def __init__(self,file_path=None):
		self.file_path = file_path

		if self.file_path != None:
			if not os.path.exists(self.file_path):
				sys.stderr.write('file {} does not exist.'.format(self.file_path))
				quit()
			try: 
				self.hdu = fits.open(self.file_path)
			except OSError as e:
				raise 

	def read_setup_yamlfile(self,setup_yamlfile):
		self.setup_yamlfile = setup_yamlfile

		if not os.path.exists(setup_yamlfile):
			raise FileNotFoundError("{} not found.".format(setup_yamlfile))
		try:
			self.param = yaml.load(open(self.setup_yamlfile))
		except OSError as e:
			raise 

	def read_xray_fitsfile(self,xray_fitsfile):
		self.xray_fitsfile = xray_fitsfile

		if not os.path.exists(self.xray_fitsfile):
			raise FileNotFoundError("{} not found.".format(self.xray_fitsfile))
		try:
			self.xray_hdu = fits.open(self.xray_fitsfile)
		except OSError as e:
			raise 
		print("{} is successfully loaded.".format(self.xray_fitsfile))

	def read_radio_fitsfile(self,mpgrp_fitsfile,ipgrp_fitsfile,radiogti_fitsfile):
		self.mpgrp_fitsfile = mpgrp_fitsfile
		self.ipgrp_fitsfile = ipgrp_fitsfile		
		self.radiogti_fitsfile = radiogti_fitsfile

		if not os.path.exists(self.mpgrp_fitsfile):
			raise FileNotFoundError("{} not found.".format(self.mpgrp_fitsfile))
		try:
			self.mpgrp_hdu = fits.open(self.mpgrp_fitsfile)
		except OSError as e:
			raise 
		print("{} is successfully loaded.".format(self.mpgrp_fitsfile))

		if not os.path.exists(self.ipgrp_fitsfile):
			raise FileNotFoundError("{} not found.".format(self.ipgrp_fitsfile))
		try:
			self.ipgrp_hdu = fits.open(self.ipgrp_fitsfile)
		except OSError as e:
			raise 
		print("{} is successfully loaded.".format(self.ipgrp_fitsfile))

		if not os.path.exists(self.radiogti_fitsfile):
			raise FileNotFoundError("{} not found.".format(self.radiogti_fitsfile))
		try:
			self.radiogti_hdu = fits.open(self.radiogti_fitsfile)
		except OSError as e:
			raise 
		print("{} is successfully loaded.".format(self.radiogti_fitsfile))

	def get_xray_profile_array(self,xray_hdu,grp_hdu,nphase=180,lag=0,flag_grp=True):
		flag_xrays_isin_grp = np.isin(
			xray_hdu['EVENTS'].data['MOD_PULSE_NUMBER']+lag,
			grp_hdu['GRP'].data['MOD_PULSE_NUMBER'])

		if flag_grp:
			print('calculation...%d (GRP)...' % lag)
			profile_count, profile_binedge, profile_patches = plt.hist(
				xray_hdu['EVENTS'].data['PULSE_PHASE'][flag_xrays_isin_grp],
				nphase,range=[0.0,1.0],histtype='step')
		else:
			print('calculation...%d (NON GRP)...' % lag)			
			profile_count, profile_binedge, profile_patches = plt.hist(
				xray_hdu['EVENTS'].data['PULSE_PHASE'][~flag_xrays_isin_grp],
				nphase,range=[0.0,1.0],histtype='step')			

		return np.array(profile_count)

	def get_xray_profile_bincenter_array(self,xray_hdu,nphase=180):
		profile_count, profile_binedge, profile_patches = plt.hist(
			xray_hdu['EVENTS'].data['PULSE_PHASE'],
			nphase,range=[0.0,1.0],histtype='step')
		profile_bincenter = 0.5*(profile_binedge[1:]+profile_binedge[:-1])
		return np.array(profile_bincenter)

	def add_header_keywords(self,hdu):
		hdu.header['DATAID'] = self.param['dataid']
		hdu.header['RADIOOBS'] = self.param['RadioObservatory']
		hdu.header['NU0'] = self.param['nu0']						
		hdu.header['NUDOT0'] = self.param['nudot0_15']*1e-15
		#hdu.header['F2'] = self.param['f2']
		hdu.header['YEAR'] = self.param['yyyy']
		hdu.header['MONTH'] = self.param['mm']				
		hdu.header['DAY'] = self.param['dd']		
		hdu.header['DOY'] = self.param['doy']				
		hdu.header['RAJ2000'] = self.param['CRAB_PULSAR_RA_J2000']				
		hdu.header['DECJ2000'] = self.param['CRAB_PULSAR_DEC_J2000']
		hdu.header['NIOBSID'] = self.param['nicerobsid']		
		hdu.header['METEPOCH'] = self.param['met_epoch']
		hdu.header['PERIOD0'] = self.param['P0ms']
		hdu.header['PDOT0'] = self.param['Pdot0']
		hdu.header['PLNTRYEP'] = self.param['PLANETARY_EPHEMERIS']
		hdu.header['SNTHR'] = self.param['SNthr']

	def copy_header_keywords(self,hdu_org,hdu_new):
		keywords = ['DATAID','RADIOOBS','NU0','NUDOT0','YEAR','MONTH','DAY','DOY',
			'RAJ2000','DECJ2000','NIOBSID','METEPOCH','PERIOD0','PDOT0','PLNTRYEP','SNTHR',
			'NXRAYEVT']
		for keyword in keywords:
			hdu_new.header[keyword] = hdu_org.header[keyword]

	def generate_fitsfile(self,setup_yamlfile,
		xray_fitsfile,mpgrp_fitsfile,ipgrp_fitsfile,radiogti_fitsfile,
		nphase=180,lagrange=10):
		self.read_setup_yamlfile(setup_yamlfile)
		self.read_xray_fitsfile(xray_fitsfile)
		self.read_radio_fitsfile(mpgrp_fitsfile,ipgrp_fitsfile,radiogti_fitsfile)

		if self.file_path == None:
			self.file_path = '%s_corr.fits' % os.path.splitext(os.path.basename(self.xray_fitsfile))[0]
		if os.path.dirname(self.file_path) != '' and not os.path.exists(os.path.dirname(self.file_path)):
			try:
				os.makedirs(os.path.dirname(self.file_path))
			except OSError as err:
				if err.errno!=17:
					raise
		profile_bincenter = self.get_xray_profile_bincenter_array(self.xray_hdu,nphase=nphase)

		cols_MPGRP = []
		cols_MPGRP.append(fits.Column(name='PULSE_PHASE',format='D',array=profile_bincenter))
		cols_NONMPGRP = []
		cols_NONMPGRP.append(fits.Column(name='PULSE_PHASE',format='D',array=profile_bincenter))
		sys.stdout.write('... making MPGRP Extension ...\n')
		for lag in range(-lagrange,lagrange+1):
			colname = 'COUNT_LAG{:0=+6}'.format(lag)
			cols_MPGRP.append(fits.Column(name=colname,format='K',unit='count',
				array=self.get_xray_profile_array(self.xray_hdu,self.mpgrp_hdu,nphase=nphase,lag=lag)
				))
			cols_NONMPGRP.append(fits.Column(name=colname,format='K',unit='count',
				array=self.get_xray_profile_array(self.xray_hdu,self.mpgrp_hdu,nphase=nphase,lag=lag,flag_grp=False)
				))						
	
		hdu_primary = fits.PrimaryHDU()
		hdu_table_mpgrp = fits.BinTableHDU.from_columns(cols_MPGRP,name='MPGRP')
		hdu_table_nonmpgrp = fits.BinTableHDU.from_columns(cols_NONMPGRP,name='NONMPGRP')		
		for hdu in [hdu_table_mpgrp,hdu_table_nonmpgrp]:
			hdu.header['NPHASE'] = nphase 
			hdu.header['PHSBINSZ'] = 1.0/float(nphase)
			self.add_header_keywords(hdu)
		hdu_table_mpgrp.header['NXRAYEVT'] = sum(hdu_table_mpgrp.data['COUNT_LAG+00000'])
		hdu_table_nonmpgrp.header['NXRAYEVT'] = sum(hdu_table_nonmpgrp.data['COUNT_LAG+00000'])	
		hdu_table_mpgrp.header['NMPGRP'] = len(self.mpgrp_hdu['GRP'].data)
		hdu_table_mpgrp.header['NIPGRP'] = 0 
		hdu_table_mpgrp.header['EXPXRAY'] = self.xray_hdu['EVENTS'].header['EXPOSURE']
		hdulist = fits.HDUList([hdu_primary,hdu_table_mpgrp,hdu_table_nonmpgrp])
		hdulist.writeto(self.file_path,overwrite=True)		

		return self.file_path 

	def get_normalized_profile(self,extension,lag):
		colname = 'COUNT_LAG{:0=+6}'.format(lag)
		sum_of_events = sum(self.hdu[extension].data[colname])

		normalized_count = self.hdu[extension].data[colname].astype(np.float32)/float(sum_of_events)
		normalized_error = np.sqrt(self.hdu[extension].data[colname].astype(np.float32))/float(sum_of_events)
		return normalized_count, normalized_error

	def write_fitsfile_with_normalized_extensions(self,outfitsfile=None):
		if outfitsfile == None:
			outfitsfile = self.file_path

		nphase = len(self.hdu['MPGRP'].data['PULSE_PHASE'])

		cols_MPGRP_norm = []
		cols_MPGRP_norm.append(fits.Column(name='PULSE_PHASE',format='D',array=self.hdu['MPGRP'].data['PULSE_PHASE']))
		for column in self.hdu['MPGRP'].columns:
			if column.name == 'PULSE_PHASE':
				continue 
			lag = int(column.name.replace('COUNT_LAG',''))
			normalized_count, normalized_error = self.get_normalized_profile('MPGRP',lag)

			colname = 'NORM_LAG{:0=+6}'.format(lag)
			cols_MPGRP_norm.append(fits.Column(name=colname,format='D',unit='norm_count',array=normalized_count))
			colname = 'NORM_ERR_LAG{:0=+6}'.format(lag)			
			cols_MPGRP_norm.append(fits.Column(name=colname,format='D',unit='norm_count',array=normalized_error))

		cols_NONMPGRP_norm = []
		cols_NONMPGRP_norm.append(fits.Column(name='PULSE_PHASE',format='D',array=self.hdu['NONMPGRP'].data['PULSE_PHASE']))
		for column in self.hdu['NONMPGRP'].columns:
			if column.name == 'PULSE_PHASE':
				continue 
			lag = int(column.name.replace('COUNT_LAG',''))
			normalized_count, normalized_error = self.get_normalized_profile('NONMPGRP',lag)

			colname = 'NORM_LAG{:0=+6}'.format(lag)
			cols_NONMPGRP_norm.append(fits.Column(name=colname,format='D',unit='norm_count',array=normalized_count))
			colname = 'NORM_ERR_LAG{:0=+6}'.format(lag)			
			cols_NONMPGRP_norm.append(fits.Column(name=colname,format='D',unit='norm_count',array=normalized_error))

		hdu_table_mpgrp_norm = fits.BinTableHDU.from_columns(cols_MPGRP_norm,name='MPGRP_NORM')
		hdu_table_nonmpgrp_norm = fits.BinTableHDU.from_columns(cols_NONMPGRP_norm,name='NONMPGRP_NORM')		

		for hdu in [hdu_table_mpgrp_norm,hdu_table_nonmpgrp_norm]:
			hdu.header['NPHASE'] = nphase 
			hdu.header['PHSBINSZ'] = 1.0/float(nphase)
		self.copy_header_keywords(self.hdu['MPGRP'],hdu_table_mpgrp_norm)
		self.copy_header_keywords(self.hdu['NONMPGRP'],hdu_table_nonmpgrp_norm)		

		hdu_primary = fits.PrimaryHDU()
		hdulist = fits.HDUList([hdu_primary,self.hdu['MPGRP'],self.hdu['NONMPGRP'],
			hdu_table_mpgrp_norm,hdu_table_nonmpgrp_norm])
		hdulist.writeto(outfitsfile,overwrite=True)		


	def plot_compared_pulseprofiles(self,lag=0,ymin=None,ymax=None,xmin=0.0,xmax=2.0):
		y_colname = 'NORM_LAG{:0=+6}'.format(lag)
		yerr_colname = 'NORM_ERR_LAG{:0=+6}'.format(lag)			

		if ymin == None:
			ymin = 0.95 * min(self.hdu['MPGRP_NORM'].data[y_colname])
		if ymax == None:
			ymax = 1.05 * max(self.hdu['MPGRP_NORM'].data[y_colname])

		peak_grp = max(self.hdu['MPGRP_NORM'].data[y_colname])
		peak_nongrp = max(self.hdu['NONMPGRP_NORM'].data[y_colname])		
		peak_err_grp = max(self.hdu['MPGRP_NORM'].data[yerr_colname])
		peak_err_nongrp = max(self.hdu['NONMPGRP_NORM'].data[yerr_colname])				

		enhance = peak_grp / peak_nongrp 
		significance = (peak_grp - peak_nongrp)/peak_err_grp

		title = '%d-%02d-%02d (%s) NICER:%s %s:SN%s LagX=%d' % (
			self.hdu['MPGRP_NORM'].header['YEAR'],
			self.hdu['MPGRP_NORM'].header['MONTH'],
			self.hdu['MPGRP_NORM'].header['DAY'],			
			self.hdu['MPGRP_NORM'].header['DATAID'],
			self.hdu['MPGRP_NORM'].header['NIOBSID'],
			self.hdu['MPGRP_NORM'].header['RADIOOBS'],
			self.hdu['MPGRP_NORM'].header['SNTHR'],
			lag 
			)
		legend_title = 'enhance=%.1f%% (%.1f-sigma)\nX-ray(%.1f ks):%d/%d\nMPGRP:%d IPGRP:%d' % (
			(enhance-1.0)*100.0,significance,
			float(self.hdu['MPGRP'].header['EXPXRAY'])*1e-3,			
			self.hdu['MPGRP'].header['NXRAYEVT'],
			self.hdu['NONMPGRP'].header['NXRAYEVT'],			
			self.hdu['MPGRP'].header['NMPGRP'],
			self.hdu['MPGRP'].header['NIPGRP']
			)

		plt.clf()
		fig, axes = plt.subplots(1,1,figsize=(9.6,5.6))
		plt.errorbar(
			np.hstack((self.hdu['NONMPGRP_NORM'].data['PULSE_PHASE'],self.hdu['NONMPGRP_NORM'].data['PULSE_PHASE']+1.0)),
			np.tile(self.hdu['NONMPGRP_NORM'].data[y_colname],2),
			yerr=np.tile(self.hdu['NONMPGRP_NORM'].data[yerr_colname],2),
			marker='',color='k',drawstyle='steps-mid',linewidth=1.0,
			label='No-GRP (%.3e +/- %.3e)' % (peak_nongrp,peak_err_nongrp))
		plt.errorbar(
			np.hstack((self.hdu['MPGRP_NORM'].data['PULSE_PHASE'],self.hdu['MPGRP_NORM'].data['PULSE_PHASE']+1.0)),
			np.tile(self.hdu['MPGRP_NORM'].data[y_colname],2),
			yerr=np.tile(self.hdu['MPGRP_NORM'].data[yerr_colname],2),
			marker='',color='r',drawstyle='steps-mid',linewidth=1.0,
			label='MP-GRP (%.3e +/- %.3e)' % (peak_grp,peak_err_grp))
		plt.vlines([1.0],0.0,ymax,'k',linestyles='dashed',linewidth=1.0)  
		plt.title(title)
		legend = axes.legend(loc='upper right',shadow=False,fontsize=11.0,title=legend_title)
		axes.set_xlim(xmin,xmax)
		axes.set_ylim(ymin,ymax)
		axes.set_xlabel('Pulse Phase')			
		axes.set_ylabel('Normalized count rate')
		axes.ticklabel_format(style="sci",axis="y",scilimits=(0,0))
		plt.subplots_adjust(wspace=0, hspace=0)
		plt.savefig('hoge.pdf',dpi=300)

		


