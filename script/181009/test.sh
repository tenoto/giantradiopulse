#!/bin/sh -f 

fcalc \
	ni1013010104_0mpu7_cl_crab_30sec_nibary_met_phase.evt \
	test_phase.evt \
	MOD_PULSE_NUMBER \
	"PULSE_PHASE < 0.5? PULSE_NUMBER: PULSE_NUMBER+1"