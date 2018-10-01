#!/bin/sh -f

export indir='/Users/enoto/work/soft/py2.7.14v180907/hoppy/hoppy/tests/timing/test_crab_pulsar_nicer/data'

make_crab_ephemeris_JPN.py data/radio/original/v181001/2017221/crab_JPNradio_ephemeris_MJD57974.yaml
mv crab_JPNradio_ephemeris_MJD57974.par data/xray/v181001/2017221	

photonphase --absphase --barytime --ephem DE430 \
        --orb $indir/ni1013010104.orb  \
        --outfile data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_photonphase.evt  \
        --plotfile data/xray/v181001/2017221/ni1013010104_0mpu7_cl_crab_30sec_photonphase.png \
        $indir/ni1013010104_0mpu7_cl_crab_30sec.evt \
        data/xray/v181001/2017221/crab_JPNradio_ephemeris_MJD57974.par	