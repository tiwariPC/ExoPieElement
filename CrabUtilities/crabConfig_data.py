""" In CRAB3 the configuration file is in Python language. It consists of creating a Configuration object imported from the WMCore library: """
from WMCore.Configuration import Configuration
config = Configuration()

"""  Once the Configuration object is created, it is possible to add new sections into it with corresponding parameters."""

config.section_("General")
config.General.requestName = 'dihiggs'
config.General.workArea = 'crab_20161226'
config.General.transferOutputs = True
config.General.transferLogs = True

DATAJEC='Summer16_23Sep2016PERIODV3_DATA'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'MVA-And-PFUnCorrectedMET.py'
config.JobType.pyCfgParams = ['runOnMC=False','period=PERIOD']
config.JobType.inputFiles = ['effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt','effAreasMuons_cone03_Spring15_25ns.txt',
'../TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml',
DATAJEC+'_Uncertainty_AK8PFchs.txt',
DATAJEC+'_Uncertainty_AK8PFPuppi.txt',
DATAJEC+'_Uncertainty_AK4PFchs.txt',
DATAJEC+'_L3Absolute_AK8PFchs.txt',
DATAJEC+'_L3Absolute_AK8PFPuppi.txt',
DATAJEC+'_L3Absolute_AK4PFchs.txt',
DATAJEC+'_L2Relative_AK8PFchs.txt',
DATAJEC+'_L2Relative_AK8PFPuppi.txt',
DATAJEC+'_L2Relative_AK4PFchs.txt',
DATAJEC+'_L1RC_AK8PFchs.txt',
DATAJEC+'_L1RC_AK4PFchs.txt',
DATAJEC+'_L1FastJet_AK8PFchs.txt',
DATAJEC+'_L1FastJet_AK8PFPuppi.txt',
DATAJEC+'_L1FastJet_AK4PFchs.txt',
DATAJEC+'_L2L3Residual_AK8PFchs.txt',
DATAJEC+'_L2L3Residual_AK8PFPuppi.txt',
DATAJEC+'_L2L3Residual_AK4PFchs.txt']
config.JobType.sendExternalFolder      = True

config.section_("Data")
config.Data.inputDataset = '/SingleMuon/Run2016B-PromptReco-v2/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20

#config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.lumiMask = 'Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.ignoreLocality = True


config.JobType.allowUndistributedCMSSW=True


#maxtarballsize = 50 
config.section_("Site")
#config.Site.storageSite = "T3_TW_NCU"
#config.Site.storageSite = "T2_CH_CERN"
##config.Site.storageSite = "T2_US_Wisconsin"                                                         
config.Site.storageSite = "T2_TW_NCHC"                                                                
