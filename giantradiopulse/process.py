__author__  = 'Teruaki Enoto'
__version__ = '0.01'
__date__    = '2018 October 24'
"""
HISTORY
2018-10-24 created by T.Enoto 
"""

import os 
import sys 
import shutil
import glob 
import yaml
import pandas as pd 

import giantradiopulse.radio
import giantradiopulse.xrayprofile 

class ObservationUnit():
	def __init__(self,row,param,outdir):
		print("-- ObservationUnit {} is generated.".format(row['dataid']))
		self.param = param
		for keyword in row.index:
			self.param[keyword] = row[keyword]

		self.outdir = '%s/%s' % (outdir,self.param['dataid'])
		self.param['suboutdir'] = self.outdir
		if not os.path.exists(self.param['suboutdir']):
			os.makedirs(self.outdir)

	def show_parameters(self):
		print(self.param)

	def set_gti_file(self,input_dir):
		search_filename = '%s/%s_*_orgFormat_GTI.txt' % (input_dir,self.param['dataid'])
		if len(glob.glob(search_filename)) == 0:
			sys.stderr.write('file {} not found.\n'.format(search_filename))
			return -1
		elif len(glob.glob(search_filename)) > 1:
			sys.stderr.write('file {} has multiple candidates.\n'.format(search_filename))
			return -1
		else:
			self.radio_gtifile_txt = glob.glob(search_filename)[0]
			sys.stdout.write('gti file {} is set.\n'.format(self.radio_gtifile_txt))
			return 0 

	def set_grplist_file(self,input_dir,grptype='MP'):
		search_filename = '%s/%s_*_%sGRPlistDE???SNge*_updated.txt' % (input_dir,self.param['dataid'],grptype)
		if len(glob.glob(search_filename)) == 0:
			sys.stderr.write('file {} not found.\n'.format(search_filename))
			return -1
		elif len(glob.glob(search_filename)) > 1:
			sys.stderr.write('file {} has multiple candidates.\n'.format(search_filename))
			return -1
		else:
			if grptype == 'MP':
				self.mpgrplist_txt = glob.glob(search_filename)[0]
				sys.stdout.write('mpgrplist file {} is set.\n'.format(self.mpgrplist_txt))
			elif grptype == 'IP':
				self.ipgrplist_txt = glob.glob(search_filename)[0]
				sys.stdout.write('ipgrplist file {} is set.\n'.format(self.ipgrplist_txt))				
			else:
				sys.stderr.write('invalid grptype [should be MP or IP]\n')				
				return -1 
			return 0 		

	def convert_radiogti_txt2fits(self):
		self.radio_gti = giantradiopulse.radio.GootTimeIntervalTextFile(self.radio_gtifile_txt)
		self.param['radio_gti_fitsfile'] = self.radio_gti.writeAsFitsFormat()
		shutil.move(self.param['radio_gti_fitsfile'],self.outdir)
		self.param['radio_gti_fitsfile'] = '%s/%s' % (self.outdir,self.param['radio_gti_fitsfile'])

	def convert_radiompgrp_txt2fits(self):
		self.radio_mpgrp = giantradiopulse.radio.GiantRadioPulseTextFile(self.mpgrplist_txt)
		self.param['radio_mpgrp_fitsfile'] = self.radio_mpgrp.writeAsFitsFormat()
		shutil.move(self.param['radio_mpgrp_fitsfile'],self.outdir)
		self.param['radio_mpgrp_fitsfile'] = '%s/%s' % (self.outdir,self.param['radio_mpgrp_fitsfile'])

	def convert_radioipgrp_txt2fits(self):
		self.radio_ipgrp = giantradiopulse.radio.GiantRadioPulseTextFile(self.ipgrplist_txt)
		self.param['radio_ipgrp_fitsfile'] = self.radio_ipgrp.writeAsFitsFormat()
		shutil.move(self.param['radio_ipgrp_fitsfile'],self.outdir)		
		self.param['radio_ipgrp_fitsfile'] = '%s/%s' % (self.outdir,self.param['radio_ipgrp_fitsfile'])

	def set_met_epoch(self):
		cmd  = 'rm -rf tmp_nitimeconv.out;'
		cmd += 'nitimeconv.py %d -f mjd -s tt > tmp_nitimeconv.out' % self.param['MJD']
		print(cmd);os.system(cmd)
		with open('tmp_nitimeconv.out') as f:
			for line in f:
				cols = line.split()
				if len(cols) == 0:
					continue 
				if cols[0] == 'NICER':
					met_epoch_day = float(cols[7])
		self.param['met_epoch'] = met_epoch_day + float(self.param['tJPLms'])*1e-3
		sys.stdout.write('MET epoch: {:.16f}\n'.format(self.param['met_epoch']))
		cmd  = 'rm -rf tmp_nitimeconv.out;'
		print(cmd);os.system(cmd)

	def set_nicer_data(self):
		candidates = glob.glob('%s/*/%s' % (self.param['PATH_TO_XRAY_DATADIR'],self.param['nicerobsid']))
		if len(candidates) == 0:
			sys.stdout.write('no nicer obsid\n')
			quit()
		self.param['path_to_nicer_obsid'] = candidates[0]

	def prepare_barycen_and_addphase_script(self):
		self.param['nibaryevt']  = '%s/ni%s_0mpu7_cl_nibary.evt' % (self.param['suboutdir'],self.param['nicerobsid'])
		self.param['nibarylog']  = '%s/ni%s_0mpu7_cl_nibary.log' % (self.param['suboutdir'],self.param['nicerobsid'])		
		self.param['niphaseevt'] = '%s/ni%s_0mpu7_cl_nibary_phase.evt'% (self.param['suboutdir'], self.param['nicerobsid'])
		self.param['niprofile'] = '%s/ni%s_0mpu7_cl_nibary_phase_profile.fht' % (self.param['suboutdir'],self.param['nicerobsid'])
		self.param['niproc_script'] = '%s/nibary_addphase_%s.sh' % (self.param['suboutdir'],self.param['nicerobsid'])

		# This calculation of F2 comes from the C code at the end of the
		# explanatory notes for the Jodrell ephemeris
		# f2 = 2.0*p1*p1/(p0*p0*p0)	
		p0 = float(self.param['P0ms'])*1e-3
		p1 = float(self.param['Pdot0'])
		self.param['f2'] = 2.0*p1*p1/(p0*p0*p0)	
		with open(self.param['niproc_script'],'w') as f:
			dump = """#!/bin/sh -f
barycorr \
infile=%s/xti/event_cl/ni%s_0mpu7_cl.evt.gz \
outfile=%s \
ra=%.6f dec=%.6f orbitfiles=%s/auxil/ni%s.orb.gz \
refframe=ICRS ephem=JPLEPH.%s \
>& %s
""" % (
	self.param['path_to_nicer_obsid'],self.param['nicerobsid'],
	self.param['nibaryevt'],
	self.param['CRAB_PULSAR_RA_J2000'],self.param['CRAB_PULSAR_DEC_J2000'],
	self.param['path_to_nicer_obsid'],self.param['nicerobsid'],
	self.param['PLANETARY_EPHEMERIS'].replace('DE',''),
	self.param['nibarylog'])
			f.write(dump)

			dump = """
rm -f %s
faddphase_nu.py %s %.7f %.13f \
--nudot=%6.5fe-15 --nu2dot=%.6e 
--outfits %s
""" % (self.param['niphaseevt'],self.param['nibaryevt'],
	self.param['met_epoch'],
	float(self.param['nu0']),float(self.param['nudot0_15']),self.param['f2'],
	self.param['niphaseevt'])
			f.write(dump)

			dump = """
rm -f %s
fplot_pulseprofile.py \
%s \
--outfits %s \
--nbin %d \
--colname PULSE_PHASE \
--title "Crab pulse profile ni%s all (MET)"    
""" % (self.param['niprofile'],
	self.param['niphaseevt'],
	self.param['niprofile'],
	self.param['NPHASE'],self.param['nicerobsid'])
			f.write(dump)

		cmd = 'chmod +x %s' % self.param['niproc_script']
		print(cmd);os.system(cmd)

	def run_barycen_and_addphase_script(self):
		cmd = '%s' % self.param['niproc_script']
		print(cmd);os.system(cmd)

	def write_parameter_yamlfile(self):
		outyamlfile = '%s/%s_setup.yaml' % (self.param['suboutdir'],self.param['dataid'])
		f = open(outyamlfile,'w')
		f.write(yaml.dump(self.param, default_flow_style=False))
		f.close()

	def reload_parameter_yamlfile(self,inputyamlfile):
		self.inputyamlfile = inputyamlfile

		if not os.path.exists(self.inputyamlfile):
			raise FileNotFoundError("{} not found.".format(self.inputyamlfile))
		try:
			self.param = yaml.load(open(self.inputyamlfile))
		except OSError as e:
			raise 
		print("setup yaml file {} is successfully loaded.".format(self.inputyamlfile))

	def prepare_datafiles(self):
		self.set_gti_file(self.param['PATH_TO_RADIO_DATADIR'])
		self.set_grplist_file(self.param['PATH_TO_RADIO_DATADIR'],grptype='MP')
		self.set_grplist_file(self.param['PATH_TO_RADIO_DATADIR'],grptype='IP')
		self.convert_radiogti_txt2fits()			
		self.convert_radiompgrp_txt2fits()
		self.convert_radioipgrp_txt2fits()
		self.set_met_epoch()
		self.set_nicer_data()
		self.prepare_barycen_and_addphase_script()
		self.run_barycen_and_addphase_script()		
		self.show_parameters()
		self.write_parameter_yamlfile()

	def generate_correlation_xrayprofile_fitsfile(self,nphase=60,lagrange=2,plot_ymin=None,plot_ymax=None):
		self.show_parameters()
		profile = giantradiopulse.xrayprofile.XrayProfileFitsfile()
		self.xrayprofile_fitsfile = profile.generate_correlation_fitsfile(
			self.inputyamlfile,
			self.param['niphaseevt'],
			self.param['radio_mpgrp_fitsfile'],
			self.param['radio_ipgrp_fitsfile'],
			self.param['radio_gti_fitsfile'],
			nphase=nphase,lagrange=lagrange)
		profile = giantradiopulse.xrayprofile.XrayProfileFitsfile(self.xrayprofile_fitsfile)
		profile.write_fitsfile_with_normalized_extensions()

		lag = 0
		lagstr = 'lag{:0=+6}'.format(lag)			
		profile = giantradiopulse.xrayprofile.XrayProfileFitsfile(self.xrayprofile_fitsfile)
		outpdf = '%s_%s.pdf' % (self.xrayprofile_fitsfile.replace('.fits',''),lagstr)
		profile.plot_compared_pulseprofiles(outpdf,lag=lag,xmin=0.90,xmax=1.10,ymin=plot_ymin,ymax=plot_ymax)		
		outpdf = '%s_%s_zoom.pdf' % (self.xrayprofile_fitsfile.replace('.fits',''),lagstr)
		profile.plot_compared_pulseprofiles(outpdf,lag=lag,xmin=0.00,xmax=2.00,ymin=plot_ymin,ymax=plot_ymax)

class ProcessManager():
	""" 
	:param file_path: path to a file to setup yaml file.
	"""

	def __init__(self,file_path,outdir='out/crabgrp'):
		self.file_path = file_path

		if not os.path.exists(self.file_path):
			raise FileNotFoundError("setup yaml file {} not found.".format(self.file_path))
		try:
			self.df = pd.read_csv(self.file_path,header=1,delim_whitespace=True)
		except OSError as e:
			raise 

		print("setup yaml file {} is successfully loaded.".format(self.file_path))

		self.param = yaml.load(open(self.file_path))
		self.outdir = outdir 

		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir)
		self.read_ephemeris_file()

	def read_ephemeris_file(self):
		self.df = pd.read_csv(self.param['CRAB_PULSAR_EPHEMERIS_FILE'],
			delim_whitespace=True,header=0)

		self.observationunit_list = []
		for index, row in self.df.iterrows():
			if row['proc_flag'] != True:
				sys.stdout.write('-- ObservationUnit {} is skipped.\n'.format(row['dataid']))
				continue 
			self.observationunit_list.append(ObservationUnit(row,self.param,self.outdir))			

	def prepare_datafiles(self):
		for obs in self.observationunit_list:
			obs.prepare_datafiles()

	def run_correlation_study(self,nphase=60,lagrange=2):
		for obs in self.observationunit_list:			
			setup_yaml = '%s/%s/%s_setup.yaml' % (self.outdir,obs.param['dataid'],obs.param['dataid'])
			obs.reload_parameter_yamlfile(setup_yaml)
			obs.generate_correlation_xrayprofile_fitsfile(
				nphase=self.param['NPHASE'],
				lagrange=self.param['LAGRANGE'],
				plot_ymin=self.param['CRAB_PROFILE_NORM_YMIN'],
				plot_ymax=self.param['CRAB_PROFILE_NORM_YMAX']
				)
			exit()

