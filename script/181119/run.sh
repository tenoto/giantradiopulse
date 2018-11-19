#!/bin/sh -f

# 2018-11-19
# rm -rf out/crabgrp/v181119
# giantradiopulse/cli/process.py prepare-datafiles script/181119/setup.yaml --outdir 'out/crabgrp/v181119'
giantradiopulse/cli/process.py add-grp-flag-to-xrayevents script/181119/setup.yaml out/crabgrp/v181119 --outdir 'out/crabgrp/v181119_grp'

#rm -rf out/crabgrp/v181105
#giantradiopulse/cli/process.py prepare-datafiles script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py run-correlation-study script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py add-observations script/181025/setup.yaml --outdir 'out/crabgrp/v181105'


#giantradiopulse/cli/process.py generate-srcspec script/181025/setup.yaml --outdir 'out/crabgrp/v181026'

#rm -rf out/crabgrp/v181105
#giantradiopulse/cli/process.py prepare-datafiles script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py run-correlation-study script/181025/setup.yaml --outdir 'out/crabgrp/v181026'
#giantradiopulse/cli/process.py add-observations script/181025/setup.yaml --outdir 'out/crabgrp/v181026'

#giantradiopulse/cli/process.py generate-bgdspec script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py generate-srcspec script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py generate-lightcurves script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
