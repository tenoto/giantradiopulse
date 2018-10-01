#!/bin/sh -f

export indir='/Users/enoto/work/soft/py2.7.14v180907/hoppy/hoppy/tests/timing/test_crab_pulsar_nicer/data'

rm -rf data/xray/v181001/2017221; mkdir -p data/xray/v181001/2017221

nibarytime.py \
	$indir/ni1013010104_0mpu7_cl_crab_30sec.evt \
	83.633218 22.014464 $indir/ni1013010104.orb \
	--refframe ICRS --ephem JPLEPH.430 \
	--outfits data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_bary.evt 

rm -f data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_phase.evt
faddphase_nu.py \
	data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_bary.evt  \
	57974.000000063641204 \
	29.6396012136058 \
	--nudot=-3.6870686500E-10 \
	--nu2dot=9.17318362810658e-21 \
	--outfits data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_mjd_phase.evt  \
	--flag_mjd  

rm -f data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_mjd_phase_profile.fht
fplot_pulseprofile.py \
	data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_mjd_phase.evt \
	--outfits data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_mjd_phase_profile.fht \
	--nbin 100 \
	--colname PULSE_PHASE \
	--title "Crab pulse profile (MJD)"    

# 57974 --> 113702332.816
# 57974.000000063641204 --> 113702332.816 + 0.0054986 = 113702332.8214986
rm -f data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase.evt 
faddphase_nu.py \
	data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_bary.evt  \
	113702332.8214986 \
	29.6396012136058 \
	--nudot=-3.6870686500E-10 \
	--nu2dot=9.17318362810658e-21 \
	--outfits data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase.evt  

rm -f data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase_profile.fht 
fplot_pulseprofile.py \
	data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase.evt \
	--outfits data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase_profile.fht \
	--nbin 100 \
	--colname PULSE_PHASE \
	--title "Crab pulse profile (MET)"    		

make_crab_ephemeris_JPN.py data/radio/original/v181001/2017221/crab_JPNradio_ephemeris_MJD57974.yaml
mv crab_JPNradio_ephemeris_MJD57974.par data/xray/v181001/2017221	

#make_crab_ephemeris_JPN.py data/radio/original/v181001/2017221/crab_JPNradio_ephemeris_MJD57974.yaml
#mv crab_JPNradio_ephemeris_MJD57974.par data/xray/v181001/2017221	