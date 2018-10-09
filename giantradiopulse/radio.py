__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 October 1'
"""
HISTORY
2018-10-09 modified by T.Enoto
2018-10-01 created by T.Enoto 
"""

import os 
import pandas as pd 
import astropy.io.fits as fits

class GootTimeInterval(object):
	def __init__(self):
		self.ngits = 0
		self.file_path = None

class GootTimeIntervalTextFile(GootTimeInterval):
	""" Represents GootTimeInterval in the Text format defined by Terasawa
	:param file_path: path to a file to be opened.
	"""

	def __init__(self,file_path):
		self.file_path = file_path

		if not os.path.exists(self.file_path):
			raise FileNotFoundError("{} not found.".format(self.file_path))
		try:
			self.df = pd.read_csv(self.file_path,header=1,delim_whitespace=True)
		except OSError as e:
			raise 

		print("{} is successfully loaded.".format(self.file_path))

	def writeAsFitsFormat(self,outfitsfile=None):
		if outfitsfile == None:
			outfitsfile = '%s.fits' % os.path.splitext(os.path.basename(self.file_path))[0]
		if os.path.dirname(outfitsfile) != '' and not os.path.exists(os.path.dirname(outfitsfile)):
			try:
				os.makedirs(os.path.dirname(outfitsfile))
			except OSError as err:
				if err.errno!=17:
					raise

		self.ngtis = self.df['nGTI']

		cols = []
		cols.append(fits.Column(name='nGTI',format='I',array=self.df['nGTI']))		
		cols.append(fits.Column(name='UTCstart',format='D',array=self.df['UTCstart'],unit='hhmmss.xxx'))
		cols.append(fits.Column(name='UTCend',format='D',array=self.df['UTCend'],unit='hhmmss.xxx'))
		cols.append(fits.Column(name='sodStart',format='D',array=self.df['sodStart'],unit='sec'))	
		cols.append(fits.Column(name='sodEnd',format='D',array=self.df['sodEnd'],unit='sec'))
		cols.append(fits.Column(name='MJDobsStart',format='D',array=self.df['MJDobsStart'],unit='day'))
		cols.append(fits.Column(name='MJDobsEnd',format='D',array=self.df['MJDobsEnd'],unit='day'))	
		cols.append(fits.Column(name='TDBstart',format='D',array=self.df['TDBstart'],unit='sec'))
		cols.append(fits.Column(name='TDBend',format='D',array=self.df['TDBend'],unit='sec'))
		cols.append(fits.Column(name='MJDtdbStart',format='D',array=self.df['MJDtdbStart'],unit='day'))				
		cols.append(fits.Column(name='MJDtdbEnd',format='D',array=self.df['MJDtdbEnd'],unit='day'))
		cols.append(fits.Column(name='duration_sec',format='D',array=self.df['duration_sec'],unit='sec'))	
		cols.append(fits.Column(name='accum_sec',format='D',array=self.df['accum_sec'],unit='sec'))
		cols.append(fits.Column(name='BARY_START',format='D',array=self.df['MJDtdbStart'],unit='day'))				
		cols.append(fits.Column(name='BARY_STOP',format='D',array=self.df['MJDtdbEnd'],unit='day'))		

		hdu_primary = fits.PrimaryHDU()
		hdu_table = fits.BinTableHDU.from_columns(cols,name='GTI')
		hdulist = fits.HDUList([hdu_primary,hdu_table])
		hdulist.writeto(outfitsfile,overwrite=True)		


def open_gti(file_path):
	if ".txt" in file_path:
		return GootTimeIntervalTextFile(file_path)
	else:
		raise NotImplementedError("GoottimeInterval class for this format is not implemented.")

class GiantRadioPulse(object):
	def __init__(self):
		self.nevents = 0
		self.file_path = None

class GiantRadioPulseTextFile(GiantRadioPulse):
	""" Represents GiantRadioPulse in the Text format defined by Terasawa
	:param file_path: path to a file to be opened.
	"""	

	def __init__(self,file_path):
		self.file_path = file_path

		if not os.path.exists(self.file_path):
			raise FileNotFoundError("{} not found.".format(self.file_path))
		try:
			self.df = pd.read_csv(self.file_path,header=4,delim_whitespace=True)
		except OSError as e:
			raise

	def writeAsFitsFormat(self,outfitsfile=None):
		if outfitsfile == None:
			outfitsfile = '%s.fits' % os.path.splitext(os.path.basename(self.file_path))[0]
		if os.path.dirname(outfitsfile) != '' and not os.path.exists(os.path.dirname(outfitsfile)):
			try:
				os.makedirs(os.path.dirname(outfitsfile))
			except OSError as err:
				if err.errno!=17:
					raise

		cols = []
		cols.append(fits.Column(name='nGRP',format='I',array=self.df['nGRP']))		
		cols.append(fits.Column(name='nORG',format='I',array=self.df['nORG']))
		cols.append(fits.Column(name='NSEQpulse',format='K',array=self.df['NSEQpulse']))
		cols.append(fits.Column(name='TDBsec',format='D',array=self.df['TDBsec'],unit='sec'))	
		cols.append(fits.Column(name='MJD',format='D',array=self.df['MJD'],unit='day'))
		cols.append(fits.Column(name='peakSN',format='D',array=self.df['peakSN']))
		cols.append(fits.Column(name='sumSN',format='D',array=self.df['sumSN']))	
		cols.append(fits.Column(name='Nsubpulse',format='I',array=self.df['Nsubpulse']))
		cols.append(fits.Column(name='phase',format='D',array=self.df['phase']))
		cols.append(fits.Column(name='UThhmmss',format='D',array=self.df['UThhmmss'],unit='yymmss.ssssss'))				
		cols.append(fits.Column(name='UTsec',format='D',array=self.df['UTsec'],unit='sec'))
		cols.append(fits.Column(name='BARY_MJD',format='D',array=self.df['MJD'],unit='day'))				


		self.nevents = len(self.df['nGRP'])

		try:
			self.df_header = pd.read_csv(self.file_path,header=0,delim_whitespace=True,nrows=1)
		except OSError as e:
			raise

		hdu_primary = fits.PrimaryHDU()
		hdu_table = fits.BinTableHDU.from_columns(cols,name='GRP')
		for colname in self.df_header.columns.values:
			if colname == 'nudot0_15':
				hdu_table.header['nudot015'] = self.df_header[colname][0]								
			else:
				hdu_table.header[colname] = self.df_header[colname][0]				
		hdulist = fits.HDUList([hdu_primary,hdu_table])
		hdulist.writeto(outfitsfile,overwrite=True)		                                                


def open_gitantradiopulse(file_path):
	if ".txt" in file_path:
		return GiantRadioPulseTextFile(file_path)
	else:
		raise NotImplementedError("GiantRadioPulse class for this format is not implemented.")

