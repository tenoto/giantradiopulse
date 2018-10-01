#!/usr/bin/env python

import os 
import giantradiopulse.radio 

def test_init():
	file_path = 'data/radio/original/v181001/2017221/2017221_UsdS_DE430NEW2_GTI.txt'
	gtifile = giantradiopulse.radio.open(file_path)

