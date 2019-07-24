import FWCore.ParameterSet.Config as cms
## removed cleaning from Exo VV package
##

process = cms.Process('NCUANA')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.options = cms.untracked.PSet(
	allowUnscheduled = cms.untracked.bool(True)
)

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.register ('runOnMC',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "runOnMC")

options.register ('isReReco',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "isReReco")

options.register ('useMiniAOD',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "useMiniAOD")

options.register ('useJECText',
		  False,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "useJECText")

options.register ('period',
		  'G',
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.string,
		  "period")


options.register ('textfiletovetoEvents',
		  'MET_Oct29/eventlist_MET_csc2015.txt',
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.string,
		  "textfiletovetoEvents")

options.parseArguments()



listEventsToSkip = []
## Only apply this for data
#if not options.runOnMC:
#    fileEventsToSkip = open(options.textfiletovetoEvents,"r")
#    for line in fileEventsToSkip:
#        cleanLine = line.rstrip()
#        listEventsToSkip.append(cleanLine+"-"+cleanLine)

#print listEventsToSkip

MCJEC='Summer16_23Sep2016V3_MC'
DATAJEC='Summer16_23Sep2016'+options.period+'V3_DATA'

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# Other statements
if options.runOnMC:
### Needs to be updated
	process.GlobalTag.globaltag='94X_mc2017_realistic_v12'
else:
    #process.GlobalTag.globaltag='92X_dataRun2_Prompt_v11'  #Conditions for prompt Prompt GT
    process.GlobalTag.globaltag='94X_dataRun2_ReReco_EOY17_v6'   #Conditions for the data reprocessing Rereco_GT
    #process.GlobalTag.globaltag='94X_dataRun2_v6'   #recommended here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD#2017_Data_re_miniAOD_31Mar2018_9



'''
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
'''

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

'''
##___________________________HCAL_Noise_Filter________________________________||
process.load('CommonTools.RecoAlgos.HBHENoiseFilterResultProducer_cfi')
process.HBHENoiseFilterResultProducer.minZeros = cms.int32(99999)
process.HBHENoiseFilterResultProducer.IgnoreTS4TS5ifJetInLowBVRegion=cms.bool(False)
#process.HBHENoiseFilterResultProducer.defaultDecision = cms.string("HBHENoiseFilterResultRun2Loose")
'''
##process.ApplyBaselineHBHENoiseFilter = cms.EDFilter('BooleanFlagFilter',
##   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHENoiseFilterResult'),
##   reverseDecision = cms.bool(False)
##)
##
##process.ApplyBaselineHBHEIsoNoiseFilter = cms.EDFilter('BooleanFlagFilter',
##   inputLabel = cms.InputTag('HBHENoiseFilterResultProducer','HBHEIsoNoiseFilterResult'),
##   reverseDecision = cms.bool(False)
##)
##


# Input source
if options.runOnMC:
	testFile='/store/mc/RunIIFall17MiniAOD/GluGluHToBB_M125_13TeV_powheg_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/20000/C8932584-5006-E811-9840-141877410512.root'#'/store/mc/RunIIFall17MiniAOD/QCD_Pt_120to170_TuneCP5_13TeV_pythia8/MINIAODSIM/94X_mc2017_realistic_v10-v1/00000/16E915A2-E60E-E811-AD53-001E67A3EF70.root'
else:
	testFile='/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/16963797-0937-E811-ABE2-008CFAE45134.root'


process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(testFile),
			    #skipEvents = cms.untracked.uint32(0)
                            )




from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsAK8'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   rParam = 0.8,
   jetCorrections = ('AK8PFchs', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
   btagDiscriminators = [
      'pfBoostedDoubleSecondaryVertexAK8BJetTags',
      'pfDeepDoubleBJetTags:probQ',
      'pfDeepDoubleBJetTags:probH',
      'pfDeepDoubleBvLJetTags:probQCD',
      'pfDeepDoubleBvLJetTags:probHbb',
      'pfDeepDoubleCvLJetTags:probQCD',
      'pfDeepDoubleCvLJetTags:probHcc',
      'pfDeepDoubleCvBJetTags:probHbb',
      'pfDeepDoubleCvBJetTags:probHcc',
      'pfMassIndependentDeepDoubleBvLJetTags:probQCD',
      'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
      'pfMassIndependentDeepDoubleCvLJetTags:probQCD',
      'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
      'pfMassIndependentDeepDoubleCvBJetTags:probHbb',
      'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
      ]
        )


## skip the events
## using the MET tails events
if options.runOnMC:
    print "No events to skip"
else:
    rangeEventsToSkip = cms.untracked.VEventRange(listEventsToSkip)
    process.source.eventsToSkip = rangeEventsToSkip

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
    process,
    isData = True, # false for MC
    fixEE2017 = True,
    fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
    postfix = "ModifiedMET"
    )

##
## This is for Uncorrected MET
from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment
## Uncorrected MET edns here
##

pvSource = 'offlineSlimmedPrimaryVertices'


bTagDiscriminators = [
    'pfJetBProbabilityBJetTags'
    ,'pfJetProbabilityBJetTags'
    ,'pfPositiveOnlyJetBProbabilityBJetTags'
    ,'pfPositiveOnlyJetProbabilityBJetTags'
    ,'pfNegativeOnlyJetBProbabilityBJetTags'
    ,'pfNegativeOnlyJetProbabilityBJetTags'
    ,'pfTrackCountingHighPurBJetTags'
    ,'pfTrackCountingHighEffBJetTags'
    ,'pfNegativeTrackCountingHighPurBJetTags'
    ,'pfNegativeTrackCountingHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighEffBJetTags'
    ,'pfSimpleSecondaryVertexHighPurBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighEffBJetTags'
    ,'pfNegativeSimpleSecondaryVertexHighPurBJetTags'
    ,'pfCombinedSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedSecondaryVertexV2BJetTags'
    ,'pfCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfPositiveCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'pfNegativeCombinedInclusiveSecondaryVertexV2BJetTags'
    ,'softPFMuonBJetTags'
    ,'positiveSoftPFMuonBJetTags'
    ,'negativeSoftPFMuonBJetTags'
    ,'softPFElectronBJetTags'
    ,'positiveSoftPFElectronBJetTags'
    ,'negativeSoftPFElectronBJetTags'
    ,'pfDeepCSVJetTags:probb'
    ,'pfDeepCSVJetTags:probc'
    ,'pfDeepCSVJetTags:probudsg'
    ,'pfDeepCSVJetTags:probbb'
]

## Jet energy corrections

## For jet energy correction
if options.runOnMC:
	jetCorrectionsAK4CHS       = ('AK4PFchs', ['L1FastJet','L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK4Puppi     = ('AK4PFPuppi', ['L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8CHS       = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8CHSL23    = ('AK8PFchs', ['L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8Puppi     = ('AK8PFPuppi', ['L2Relative', 'L3Absolute'], 'None')
	jetCorrectionLevelsFullCHS = ['L1FastJet', 'L2Relative', 'L3Absolute']
	jetCorrectionLevels23CHS   = ['L2Relative', 'L3Absolute']
	jetCorrectionLevelsPuppi   = ['L2Relative', 'L3Absolute']

	AK4JECTextFiles = [
		MCJEC+'_L1FastJet_AK4PFchs.txt',
		MCJEC+'_L2Relative_AK4PFchs.txt',
		MCJEC+'_L3Absolute_AK4PFchs.txt'
		]
	AK4JECUncTextFile = MCJEC+'_Uncertainty_AK4PFchs.txt'

	AK8JECTextFiles = [
		MCJEC+'_L1FastJet_AK8PFchs.txt',
		MCJEC+'_L2Relative_AK8PFchs.txt',
		MCJEC+'_L3Absolute_AK8PFchs.txt'
		]
	AK8JECUncTextFile = MCJEC+'_Uncertainty_AK8PFchs.txt'
	prunedMassJECTextFiles = [
		MCJEC+'_L2Relative_AK8PFchs.txt',
		MCJEC+'_L3Absolute_AK8PFchs.txt'
		]

	AK4PuppiJECTextFiles = [
		MCJEC+'_L2Relative_AK4PFPuppi.txt',
		MCJEC+'_L3Absolute_AK4PFPuppi.txt'
		]
	AK4PuppiJECUncTextFile = MCJEC+'_Uncertainty_AK4PFPuppi.txt'

	AK8PuppiJECTextFiles = [
		MCJEC+'_L2Relative_AK8PFPuppi.txt',
		MCJEC+'_L3Absolute_AK8PFPuppi.txt'
		]
	AK8PuppiJECUncTextFile = MCJEC+'_Uncertainty_AK8PFPuppi.txt'
else:
        jetCorrectionsAK4CHS       = ('AK4PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK4Puppi     = ('AK4PFPuppi', ['L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8CHS       = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8CHSL23    = ('AK8PFchs', ['L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8Puppi     = ('AK8PFPuppi', ['L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionLevelsFullCHS = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
	jetCorrectionLevels23CHS   = ['L2Relative', 'L3Absolute','L2L3Residual']
	jetCorrectionLevelsPuppi   = ['L2Relative', 'L3Absolute','L2L3Residual']
	AK4JECTextFiles = [
		DATAJEC+'_L1FastJet_AK4PFchs.txt',
		DATAJEC+'_L2Relative_AK4PFchs.txt',
		DATAJEC+'_L3Absolute_AK4PFchs.txt',
		DATAJEC+'_L2L3Residual_AK4PFchs.txt'
		]
	AK4JECUncTextFile = DATAJEC+'_Uncertainty_AK4PFchs.txt'
	AK8JECTextFiles = [
		DATAJEC+'_L1FastJet_AK8PFchs.txt',
		DATAJEC+'_L2Relative_AK8PFchs.txt',
		DATAJEC+'_L3Absolute_AK8PFchs.txt',
		DATAJEC+'_L2L3Residual_AK8PFchs.txt'
		]
	AK8JECUncTextFile = DATAJEC+'_Uncertainty_AK8PFchs.txt'
	prunedMassJECTextFiles = [
		DATAJEC+'_L2Relative_AK8PFchs.txt',
		DATAJEC+'_L3Absolute_AK8PFchs.txt',
		DATAJEC+'_L2L3Residual_AK8PFchs.txt'
		]

	AK4PuppiJECTextFiles = [
		DATAJEC+'_L2Relative_AK4PFPuppi.txt',
		DATAJEC+'_L3Absolute_AK4PFPuppi.txt',
		DATAJEC+'_L2L3Residual_AK4PFPuppi.txt'
		]
	AK4PuppiJECUncTextFile = DATAJEC+'_Uncertainty_AK4PFPuppi.txt'

	AK8PuppiJECTextFiles = [
		DATAJEC+'_L2Relative_AK8PFPuppi.txt',
		DATAJEC+'_L3Absolute_AK8PFPuppi.txt',
		DATAJEC+'_L2L3Residual_AK8PFPuppi.txt'
		]
	AK8PuppiJECUncTextFile = DATAJEC+'_Uncertainty_AK8PFPuppi.txt'

from PhysicsTools.PatAlgos.tools.jetTools import *


NOTADDHBBTag=False
## Filter for good primary vertex
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
    vertexCollection = cms.InputTag(pvSource),
    minimumNDOF = cms.uint32(4) ,
    maxAbsZ = cms.double(24),
    maxd0 = cms.double(2)
)

#-------------------------------------
from PhysicsTools.PatAlgos.tools.pfTools import *
## Adapt primary vertex collection
adaptPVs(process, pvCollection=cms.InputTag(pvSource))



#### Add reclustered AK8 Puppi jet by Eiko


process.load('CommonTools/PileupAlgos/Puppi_cff')
process.puppi.candName       = cms.InputTag('packedPFCandidates')
process.puppi.vertexName     = cms.InputTag('offlineSlimmedPrimaryVertices')

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox


### CA15Puppi
jetToolbox( process, 'ca15', 'jetSequence', 'out', PUMethod='Puppi', miniAOD=options.useMiniAOD, runOnMC=options.runOnMC,
	    bTagDiscriminators=(bTagDiscriminators + ([] if NOTADDHBBTag else ['pfBoostedDoubleSecondaryVertexCA15BJetTags'])),
	    JETCorrPayload='AK8PFPuppi',JETCorrLevels=jetCorrectionLevelsPuppi,
	    subJETCorrPayload='AK4PFPuppi',subJETCorrLevels=jetCorrectionLevelsPuppi,
	    Cut='pt>120',
	    addSoftDrop=True,addSoftDropSubjets=True, betaCut=1.0, zCutSD=0.15,
	    addNsub=True )

### AK8Puppi
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='Puppi', miniAOD=options.useMiniAOD, runOnMC=options.runOnMC,
#	    bTagDiscriminators=(bTagDiscriminators + ([] if NOTADDHBBTag else ['pfBoostedDoubleSecondaryVertexAK8BJetTags'])),
#	    JETCorrPayload='AK8PFPuppi',JETCorrLevels=jetCorrectionLevelsPuppi,
#	    subJETCorrPayload='AK4PFPuppi',subJETCorrLevels=jetCorrectionLevelsPuppi,
#	    Cut='pt>170',
#	    addSoftDrop=True,addSoftDropSubjets=True,addNsub=True )


## and add them to the event content

'''
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContentNoCleaning
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('patTuple.root'),
                               ## save only events passing the full path
                               #SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               ## save PAT output; you need a '*' to unpack the list of commands
                               ## 'patEventContent'
                               outputCommands = cms.untracked.vstring('drop *', *patEventContentNoCleaning )
                               )

#process.myLabel = cms.EDAnalyzer('DemoAnalyzer',
#            JetTag      = cms.InputTag('selectedUpdatedPatJets'),
#)


patAlgosToolsTask = getPatAlgosToolsTask(process)
process.outpath = cms.EndPath(process.out, patAlgosToolsTask)

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsAK8'),
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   rParam = 0.8,
   jetCorrections = ('AK8PFchs', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
   btagDiscriminators = [
      'pfBoostedDoubleSecondaryVertexAK8BJetTags',
      'pfDeepDoubleBJetTags:probQ',
      'pfDeepDoubleBJetTags:probH',
      'pfDeepDoubleBvLJetTags:probQCD',
      'pfDeepDoubleBvLJetTags:probHbb',
      'pfDeepDoubleCvLJetTags:probQCD',
      'pfDeepDoubleCvLJetTags:probHcc',
      'pfDeepDoubleCvBJetTags:probHbb',
      'pfDeepDoubleCvBJetTags:probHcc',
      'pfMassIndependentDeepDoubleBvLJetTags:probQCD',
      'pfMassIndependentDeepDoubleBvLJetTags:probHbb',
      'pfMassIndependentDeepDoubleCvLJetTags:probQCD',
      'pfMassIndependentDeepDoubleCvLJetTags:probHcc',
      'pfMassIndependentDeepDoubleCvBJetTags:probHbb',
      'pfMassIndependentDeepDoubleCvBJetTags:probHcc',
      ]
        )

from Configuration.EventContent.EventContent_cff import MINIAODSIMEventContent
process.out.outputCommands.append('keep *_slimmedJetsAK8*_*_*')
process.out.outputCommands.append('keep *_offlineSlimmedPrimaryVertices*_*_*')
process.out.outputCommands.append('keep *_slimmedSecondaryVertices*_*_*')
process.out.outputCommands.append('keep *_selectedPatJets*_*_*')
process.out.outputCommands.append('keep *_selectedUpdatedPatJets*_*_*')
process.out.outputCommands.append('keep *_pfBoostedDoubleSVAK8TagInfos*_*_*')
process.out.outputCommands.append('keep *_pfDeepDoubleXTagInfos*_*_*')
process.out.outputCommands.append('keep *_updatedPatJets*_*_*')
'''

###end of add jet collection



## add value maps for electron IDs
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate
if options.useMiniAOD:
	dataFormat = DataFormat.MiniAOD
else :
	dataFormat = DataFormat.AOD

switchOnVIDElectronIdProducer(process, dataFormat)
switchOnVIDPhotonIdProducer(process, dataFormat)
# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V1_cff', 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V1_cff']


#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V1_cff', 'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V1_cff',]

#add them to the VID producer
for idmod in my_phoid_modules:
	setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)




#process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("ncuslimmedPhoton")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("ncuslimmedElectron")

from DelPanj.TreeMaker.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
    debug=True,
    toKeep = ["2017v1", "2017v2", "newDM2017v2", "dR0p32017v2", "2016v1", "newDM2016v1"]
)
na.runTauID()

byIsolationMVArun2017v2DBoldDMwLTraw2017 = cms.string('byIsolationMVArun2017v2DBoldDMwLTraw2017'),
byVVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVVLooseIsolationMVArun2017v2DBoldDMwLT2017'),
byVLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVLooseIsolationMVArun2017v2DBoldDMwLT2017'),
byLooseIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byLooseIsolationMVArun2017v2DBoldDMwLT2017'),
byMediumIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byMediumIsolationMVArun2017v2DBoldDMwLT2017'),
byTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byTightIsolationMVArun2017v2DBoldDMwLT2017'),
byVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVTightIsolationMVArun2017v2DBoldDMwLT2017'),
byVVTightIsolationMVArun2017v2DBoldDMwLT2017 = cms.string('byVVTightIsolationMVArun2017v2DBoldDMwLT2017')

## For normal AK4 jets jet energy correction on top of miniAOD
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJECAK4 = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("appliedRegJets"),
	levels = jetCorrectionLevelsFullCHS,
	payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJECAK4 = updatedPatJets.clone(
	jetSource = cms.InputTag("appliedRegJets"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK4"))
  )

process.jetCorrSequenceAK4 = cms.Sequence( process.patJetCorrFactorsReapplyJECAK4 + process.patJetsReapplyJECAK4 )




### For normal AK8 jet energy correction on top of miniAOD
process.patJetCorrFactorsReapplyJECAK8 = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJetsAK8"),
	levels = jetCorrectionLevelsFullCHS,
	payload = 'AK8PFPuppi' ) # Make sure to choose the appropriate levels and payload here!

process.patJetsReapplyJECAK8 = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJetsAK8"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8"))
  )


process.jetCorrSequenceAK8 = cms.Sequence( process.patJetCorrFactorsReapplyJECAK8 + process.patJetsReapplyJECAK8 )



## For normal AK4Puppi jets jet energy correction on top of miniAOD
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJECAK4Puppi = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJetsPuppi"),
	levels = jetCorrectionLevelsPuppi,
	payload = 'AK4PFPuppi' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJECAK4Puppi = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJetsPuppi"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK4Puppi"))
  )

process.jetCorrSequenceAK4Puppi = cms.Sequence( process.patJetCorrFactorsReapplyJECAK4Puppi + process.patJetsReapplyJECAK4Puppi )


## For correcting pruned jet mass + CHS
process.patJetCorrFactorsReapplyJECForPrunedMass = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJetsAK8"),
	levels = jetCorrectionLevels23CHS,
	payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

process.patJetsReapplyJECForPrunedMass = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJetsAK8"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECForPrunedMass"))
	)

process.jetCorrSequenceForPrunedMass = cms.Sequence( process.patJetCorrFactorsReapplyJECForPrunedMass + process.patJetsReapplyJECForPrunedMass )




process.load('DelPanj.TreeMaker.TreeMaker_cfi')
process.tree.useJECText            = cms.bool(options.useJECText)
process.tree.THINjecNames          = cms.vstring(AK4JECTextFiles)
process.tree.THINjecUncName        = cms.string(AK4JECUncTextFile)
process.tree.FATprunedMassJecNames = cms.vstring(prunedMassJECTextFiles)
process.tree.FATjecNames           = cms.vstring(AK8PuppiJECTextFiles)
process.tree.FATjecUncName         = cms.string(AK8PuppiJECUncTextFile)
process.tree.AK4PuppijecNames      = cms.vstring(AK4PuppiJECTextFiles)
process.tree.AK4PuppijecUncName    = cms.string(AK4PuppiJECUncTextFile)
process.tree.AK8PuppijecNames      = cms.vstring(AK8PuppiJECTextFiles)
process.tree.AK8PuppijecUncName    = cms.string(AK8PuppiJECUncTextFile)
process.tree.CA15PuppijecNames     = cms.vstring(AK8PuppiJECTextFiles)
process.tree.CA15PuppijecUncName   = cms.string(AK8PuppiJECUncTextFile)
process.tree.fillCA15PuppiJetInfo  = cms.bool(True)


if options.useJECText:
	process.tree.THINJets      = cms.InputTag("appliedRegJets")
	process.tree.FATJets       = cms.InputTag("selectedUpdatedPatJets")#("slimmedJetsAK8")
	process.tree.FATJetsForPrunedMass       = cms.InputTag("slimmedJetsAK8")
	process.tree.AK4PuppiJets  = cms.InputTag("slimmedJetsPuppi")



process.TFileService = cms.Service("TFileService",fileName = cms.string("NCUGlobalTuples.root"))

##Trigger Filter
process.trigFilter = cms.EDFilter('TrigFilter',
                                 TrigTag = cms.InputTag("TriggerResults::HLT"),
                                 TrigPaths = cms.vstring("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight","HLT_PFMETNoMu140_PFMHTNoMu140_IDTight","HLT_IsoMu27","HLT_IsoTkMu27","HLT_IsoMu24","HLT_IsoTkMu24","HLT_Ele27_WPTight_Gsf","HLT_Ele32_WPTight_Gsf_L1DoubleEG","HLT_Ele35_WPTight_Gsf","HLT_Photon200"),
                                 isMC_ = cms.bool(options.runOnMC)
                                 )
## New MET Filters
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadPFMuonFilter.taggingMode = cms.bool(True)

process.load('RecoMET.METFilters.badGlobalMuonTaggersMiniAOD_cff')
process.badGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)
process.cloneGlobalMuonTaggerMAOD.taggingMode = cms.bool(True)

##
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")
process.BadChargedCandidateFilter.taggingMode = cms.bool(True)
##
##
process.allEventsCounter = cms.EDFilter(
	"EventCounter"
 )



process.appliedRegJets= cms.EDProducer('bRegressionProducer',
                                           JetTag=cms.InputTag("slimmedJets"),
                                           rhoFixedGridCollection = cms.InputTag('fixedGridRhoFastjetAll'),
                                           #bRegressionWeightfile= cms.untracked.string("/afs/cern.ch/work/d/dekumar/public/flashgg_setup/CMSSW_8_0_28/src/flashgg/MetaData/data/DNN_models/model-18"),
                                           y_mean = cms.untracked.double(1.0454729795455933) ,
                                           y_std = cms.untracked.double( 0.31628304719924927)
                                           )

if not options.useJECText:
	process.analysis = cms.Path(
	process.appliedRegJets+  #added by deepak
        process.trigFilter+
		process.allEventsCounter+
		process.egmGsfElectronIDSequence+
		process.egmPhotonIDSequence+
		process.rerunMvaIsolationSequence
		*process.NewTauIDsEmbedded+
		process.fullPatMetSequenceModifiedMET+
		process.pfMet+
		process.jetCorrSequenceAK4+
		process.jetCorrSequenceAK8+
		process.jetCorrSequenceAK4Puppi+
		process.jetCorrSequenceForPrunedMass+
		process.BadPFMuonFilter +
		process.BadChargedCandidateFilter +
		process.badGlobalMuonTaggerMAOD +
		process.cloneGlobalMuonTaggerMAOD +
		#process.HBHENoiseFilterResultProducer+ ## by raman
		process.tree
		)
else:
	process.analysis = cms.Path(
	process.appliedRegJets+
        process.trigFilter+
		process.allEventsCounter+
		process.egmGsfElectronIDSequence+
		process.egmPhotonIDSequence+
		process.rerunMvaIsolationSequence*
		process.NewTauIDsEmbedded+
		process.fullPatMetSequenceModifiedMET+
		process.pfMet+
		process.BadPFMuonFilter+
		process.BadChargedCandidateFilter+
		process.badGlobalMuonTaggerMAOD+
		process.cloneGlobalMuonTaggerMAOD+
		process.tree
		)


#print process.dumpPython()
