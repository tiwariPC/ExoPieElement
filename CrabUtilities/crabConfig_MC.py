""" In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: """
from WMCore.Configuration import Configuration
config = Configuration()

"""  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters."""

config.section_("General")
config.General.requestName = 'dihiggs'
config.General.workArea = 'crab_20161226'
config.General.transferOutputs = True
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'MVA-And-PFUnCorrectedMET.py'
config.JobType.inputFiles = ['effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt','effAreasMuons_cone03_Spring15_25ns.txt',
'Summer16_23Sep2016V3_MC_Uncertainty_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_Uncertainty_AK8PFPuppi.txt',
'Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt',
'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_L3Absolute_AK8PFPuppi.txt',
'Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt',
'Summer16_23Sep2016V3_MC_L2Relative_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_L2Relative_AK8PFPuppi.txt',
'Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt',
'Summer16_23Sep2016V3_MC_L2L3Residual_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_L2L3Residual_AK8PFPuppi.txt',
'Summer16_23Sep2016V3_MC_L2L3Residual_AK4PFchs.txt',
'Summer16_23Sep2016V3_MC_L1RC_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_L1RC_AK4PFchs.txt',
'Summer16_23Sep2016V3_MC_L1FastJet_AK8PFchs.txt',
'Summer16_23Sep2016V3_MC_L1FastJet_AK8PFPuppi.txt',
'Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt',
'../TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml']
config.JobType.sendExternalFolder      = True

config.section_("Data")
config.Data.inputDataset = '/DYJetsToLL_M-50_HT-100to200_Tune4C_13TeV-madgraph-tauola/Phys14DR-PU20bx25_PHYS14_25_V1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
##config.Data.outLFNDirBase = '/store/user/khurana/MonoH2016/V3/'


config.JobType.allowUndistributedCMSSW=True


config.section_("Site")
config.Site.storageSite = "T3_TW_NCU"
#config.Site.storageSite = "T2_CH_CERN"
##config.Site.storageSite = "T2_US_Wisconsin"
#config.Site.storageSite = "T2_TW_NCHC"
