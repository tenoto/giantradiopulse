#!/bin/sh -f

#rm -rf out/crabgrp/v181105
#giantradiopulse/cli/process.py prepare-datafiles script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
giantradiopulse/cli/process.py run-correlation-study script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py add-observations script/181025/setup.yaml --outdir 'out/crabgrp/v181105'


#giantradiopulse/cli/process.py generate-srcspec script/181025/setup.yaml --outdir 'out/crabgrp/v181026'

#rm -rf out/crabgrp/v181105
#giantradiopulse/cli/process.py prepare-datafiles script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py run-correlation-study script/181025/setup.yaml --outdir 'out/crabgrp/v181026'
#giantradiopulse/cli/process.py add-observations script/181025/setup.yaml --outdir 'out/crabgrp/v181026'

#giantradiopulse/cli/process.py generate-bgdspec script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py generate-srcspec script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
#giantradiopulse/cli/process.py generate-lightcurves script/181025/setup.yaml --outdir 'out/crabgrp/v181105'
