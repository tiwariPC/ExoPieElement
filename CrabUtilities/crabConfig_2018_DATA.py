from WMCore.Configuration import Configuration
config = Configuration()

workname='setup_2018_v07'
reqname='DYJetsToLL_M_50_HT_400to600_TuneCP5_13TeV_30K' ## this is temporrt, will be changed by the multicrab
dataset='/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v4/MINIAODSIM' ## this is temporrt, will be changed by the multicrab
number_of_units=1 ## this is temporrt, will be changed by the multicrab

config.section_("General")
config.General.requestName = 'XYZ_replacedbyMulticrab'
config.General.workArea = 'crab_'+workname
config.General.transferOutputs = True
config.General.transferLogs = False

PERIOD="D"
DATAJEC='Summer16_23Sep2016'+PERIOD+'V3_DATA' ## this is not used at this moment

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeMaker_2018_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=False','runOn2018=True']
config.JobType.inputFiles = ['../MetaData/data/DNN_models/model-37.pb',
'../TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml',
'../TreeMaker/data/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt']

config.JobType.sendExternalFolder = True
config.JobType.sendPythonFolder = True
config.JobType.outputFiles = ['ExoPieElementTuples.root']
config.JobType.disableAutomaticOutputCollection = True

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.outputDatasetTag = reqname

#config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.splitting = 'EventAwareLumiBased'
#config.Data.unitsPerJob = 4000

#config.Data.ignoreLocality = True
config.JobType.allowUndistributedCMSSW=True


#maxtarballsize = 50
config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"

config.Data.outLFNDirBase = '/store/group/phys_exotica/bbMET/setup_2018_v07'