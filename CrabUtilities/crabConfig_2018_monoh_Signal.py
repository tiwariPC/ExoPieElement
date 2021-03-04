""" In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: """
from WMCore.Configuration import Configuration
config = Configuration()

"""  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters."""

workname='setup_2018_v05'
reqname='scalar_NLO_Mchi-450_Mphi-1000'
dataset='/BBbarDMJets_scalar_NLO_Mchi-450_Mphi-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'

PERIOD='H'
EPLS=10


config.section_("General")
config.General.requestName = reqname
config.General.workArea = 'crab_'+workname
config.General.transferOutputs = True
config.General.transferLogs = False

#DATAJEC='Summer16_23Sep2016'+PERIOD+'V3_DATA'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeMaker_2018_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=True','runOn2018=True','runOnrpZpB=True']
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
#config.Data.outputDatasetTag = reqname

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
