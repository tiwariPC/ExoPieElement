""" In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: """
from WMCore.Configuration import Configuration
config = Configuration()

"""  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters."""

workname='2017_data_20190522'
reqname='SingleElectron-Run2017C-31Mar2018-v1'
dataset='/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD'

PERIOD='H'

config.section_("General")
config.General.requestName = reqname
config.General.workArea = 'crab_'+workname
config.General.transferOutputs = True
config.General.transferLogs = True

DATAJEC='Summer16_23Sep2016'+PERIOD+'V3_DATA'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'treeMaker_Summer17_cfg.py'
config.JobType.pyCfgParams = ['runOnMC=False','period='+PERIOD]
config.JobType.inputFiles = ['effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt','effAreasMuons_cone03_Spring15_25ns.txt',
'../../../MetaData/data/DNN_models/breg_training_2017.pb',
'../../TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml']

config.JobType.sendExternalFolder = True
config.JobType.sendPythonFolder = True

config.section_("Data")
config.Data.inputDataset = dataset
config.Data.inputDBS = 'global'
config.Data.outputDatasetTag = reqname
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'

config.Data.unitsPerJob = 20

#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'#'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
#config.Data.ignoreLocality = True


config.JobType.allowUndistributedCMSSW=True


#maxtarballsize = 50 
config.section_("Site")
#config.Site.storageSite = "T3_TW_NCU"
config.Site.storageSite = "T2_CH_CERN"
##config.Site.storageSite = "T2_US_Wisconsin"                                                         
#config.Site.storageSite = "T2_TW_NCHC"                                                               

config.Data.outLFNDirBase = '/store/group/phys_exotica/bbMET/2017_ntuples_v2/%s' %(workname) 
