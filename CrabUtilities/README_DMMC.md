
## Check the status of step4 production and prepare input for DMAnaRun2

Go to your work directory for step3 and step4 (created following the instruction at https://github.com/syuvivida/DM/tree/v0.05/signalMC_production_crab/step3_step4)

Make sure your step4 jobs are 100% finished before you move on to DMAnaRun2

```
source /afs/cern.ch/cms/cmsset_default.csh (bash: source /afs/cern.ch/cms/cmsset_default.sh)
cd CMSSW_8_0_11/src
cmsenv
source /cvmfs/cms.cern.ch/crab3/crab_standalone.csh (bash: source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh)
curl https://raw.githubusercontent.com/syuvivida/DMAnaRun2/80X_dev/CrabUtilities/MultiCrab_DMAna.py -o MultiCrab_DMAna.py
voms-proxy-init --voms cms
python MultiCrab_DMAna.py prepare step4_inputdataset.txt
```
A text file "ana_inputdataset.txt" will be created.


## Software setup before you submit CRAB jobs

Note, you should follow the instruction here but do not run any MultiCrab or submit any CRAB jobs yet 
https://github.com/syuvivida/DMAnaRun2 

## Submit CRAB jobs

under the directory DelPanj/CrabUtilities, 

copy the file ana_inputdataset.txt to this directory

```
source /cvmfs/cms.cern.ch/crab3/crab_standalone.csh (bash: source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh)
curl https://raw.githubusercontent.com/syuvivida/DMAnaRun2/80X_dev/CrabUtilities/MultiCrab_DMAna.py -o MultiCrab_DMAna.py
voms-proxy-init --voms cms
python MultiCrab_DMAna.py submit
```

## Check status of CRAB jobs
```
source /cvmfs/cms.cern.ch/crab3/crab_standalone.csh (bash: source /cvmfs/cms.cern.ch/crab3/crab_standalone.sh)
voms-proxy-init --voms cms
python MultiCrab_DMAna.py status <dirName>
```