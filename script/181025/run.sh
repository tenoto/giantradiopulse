#!/bin/sh -f

#rm -rf out/crabgrp/v181026
#giantradiopulse/cli/process.py prepare-datafiles script/181025/setup.yaml --outdir 'out/crabgrp/v181026'
giantradiopulse/cli/process.py run-correlation-study script/181025/setup.yaml --outdir 'out/crabgrp/v181026'