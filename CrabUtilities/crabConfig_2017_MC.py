from WMCore.Configuration import Configuration
config = Configuration()

workname ='MC_2017miniaodV2_V1'
reqname='DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV_30K' ## this is temporrt, will be changed by the multicrab
dataset='/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM' ## this is temporrt, will be changed by the multicrab
number_of_units=1 ## this is temporrt, will be changed by the multicrab

config.section_("General")
config.General.requestName = 'XYZ_replacedbyMulticrab'
config.General.workArea = 'crab_'+workname
config.General.transferOutputs = True
config.General.transferLogs = True

PERIOD="D"
DATAJEC='Summer16_23Sep2016'+PERIOD+'V3_DATA' ## this is not used at this moment

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeMaker_16_17_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=True']

config.JobType.inputFiles = ['../MetaData/data/DNN_models/breg_training_2017.pb',  ## this is for b-jet regression
                             '../TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml'  ## this is for CA15 training, will be removed in next iteration
                             ]

#config.JobType.sendExternalFolder = True
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
#config.Data.outputDatasetTag = config.General.requestName

#config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = number_of_units

config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 30000

#config.Data.ignoreLocality = True
config.JobType.allowUndistributedCMSSW=True


#maxtarballsize = 50
config.section_("Site")
#config.Site.storageSite = "T3_TW_NCU"
config.Site.storageSite = "T2_CH_CERN"
##config.Site.storageSite = "T2_US_Wisconsin"
#config.Site.storageSite = "T2_TW_NCHC"

config.Data.outLFNDirBase = '/store/group/phys_exotica/bbMET/ExoPieElementTuples/%s' %(workname)
