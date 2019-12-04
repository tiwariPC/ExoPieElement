## this branch is already good for 2017 analysis, 
## add the instruction for 2018
## please make all changes for 2016 in the present branch such that this can be used for both years. 

## make all the new switches configurable via python file, DO NOT MAKE C++ ONLY SWITCHES ANYMORE. 

# Installation of ExoPieElement and dependencies
New README

## now the setup works both in SLC6 and centos 7, You can login to slc6 using USERNAME@lxplus6.cern.ch   and slc7 using USERNAME@lxplus.cern.ch

## it is recomended to use centos 7 

### setup CMSSW

export SCRAM_ARCH=slc7_amd64_gcc630     ## for slc 6 use export SCRAM_ARCH=slc6_amd64_gcc630 

cmsrel CMSSW_9_4_13

cd CMSSW_9_4_13/src

cmsenv

git cms-init

git cms-merge-topic cms-egamma:EgammaPostRecoTools #just adds in an extra file to have a setup function to make things easier


git cms-merge-topic cms-met:METFixEE2017_949_v2


git clone git@github.com:ExoPie/ExoPieElement.git

scram b -j 4 


##For jetToolBox

git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox

cd JMEAnalysis/JetToolbox

git checkout jetToolbox_94X_v3

cd -

scram b -j 4


cd ExoPieElement/TreeMaker

config file to run is in test dir: named treeMaker_16_17_cfg.py 

cd test 

## before cmsRun do set the proxy. 

cmsRun treeMaker_16_17_cfg.py

If the file doesn't work, instead of /tmp/khurana.... use filename, it will take some time to run via xrootd. 

you have to change the line 

fileNames = cms.untracked.vstring("file:/tmp/khurana/temp2017.root"),

by

fileNames = cms.untracked.vstring(testFile),


## for crab submission 

go the the directory: 

cd ExoPieElement/CrabUtilities

copy the latest .py file from test dir, do it via

source setup.sh


edit the crabConfig_2017_MC.py file, change following parameters as per your site of storage: 

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 30000

config.Site.storageSite = "T3_TW_NCU"                                                                                                                                                                      
config.Data.outLFNDirBase = '/store/group/phys_exotica/bbMET/ExoPieElementTuples/%s' %(workname)

you can make a list of samples in a text and use it to submit multiple jobs using MultiCrab_2017MC.py

you need to edit the .txt file name and run it using 

python MultiCrab_2017MC.py

 to check the status of all the jobs you just submitted you can add one function in the same file to do it. there is some example in 




for crab submit: 
python MultiCrab_2017MC.py --submit 

for crab status: 
python MultiCrab_2017MC.py --status --crabdir=crab_MC_2017miniaodV2_V1

for crab status with summary of all the dataset:  (ss refer to status summary)
python MultiCrab_2017MC.py --status --crabdir=crab_MC_2017miniaodV2_V1 --ss 

for crab resubmit: 
python MultiCrab_2017MC.py --resubmit --crabdir=crab_MC_2017miniaodV2_V1

for crab kill: 
python MultiCrab_2017MC.py --kill --crabdir=crab_MC_2017miniaodV2_V1




## Instruction for 2018 setup

Working on lxplus7

##SCRAM ARCH SET for CMSSW_10_2_10

export SCRAM_ARCH=slc7_amd64_gcc700

cmsrel CMSSW_10_2_17

cd CMSSW_10_2_17/src/

git cms-init

git cms-merge-topic cms-egamma:EgammaPostRecoTools

git cms-merge-topic cms-met:METFixEE2017_949_v2_backport_to_102X  ##Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETUncertaintyPrescription

git clone -b setup_2017_2016_2018 git@github.com:ExoPie/ExoPieElement.git

scram b -j 4

git clone git@github.com:cms-jet/JetToolbox.git JMEAnalysis/JetToolbox -b jetToolbox_102X_v2

scram b -j 4

cd ExoPieElement/TreeMaker/test/

##Then run the test file by

cmsRun treeMaker_2018_cfg.py



