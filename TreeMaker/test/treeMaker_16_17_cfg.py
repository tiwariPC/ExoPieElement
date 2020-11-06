import FWCore.ParameterSet.Config as cms

process = cms.Process('EXOPIE')
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


options.register ('period',
		    'G',
		    VarParsing.multiplicity.singleton,
		    VarParsing.varType.string,
		    "period")

options.register ('useJECText',
		  False,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "useJECText")

options.register ('useMiniAOD',
		    True,
		    VarParsing.multiplicity.singleton,
		    VarParsing.varType.bool,
		    "useMiniAOD")

options.register ('runOn2017',
		  False,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "runOn2017")

options.register ('runOn2016',
		  False,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "runOn2016")

options.parseArguments()


MCJEC='Summer16_23Sep2016V3_MC'
DATAJEC='Summer16_23Sep2016'+options.period+'V3_DATA'

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')


from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmVAnalysisSummaryTable
if options.runOn2017:
    if options.runOnMC:
        process.GlobalTag.globaltag='94X_mc2017_realistic_v17'
    else:
        process.GlobalTag.globaltag='94X_dataRun2_v11'
elif options.runOn2016:
    if options.runOnMC:
        process.GlobalTag.globaltag='94X_mcRun2_asymptotic_v3'
    else:
        process.GlobalTag.globaltag='94X_dataRun2_v10'


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

## New from Egamma
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPostRecoRecipes#Running_on_2017_MiniAOD_V2
## for 2017 recomended ID is Fall17V2
## a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True, #if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2017-Nov17ReReco')  #era is new to select between 2016 / 2017,  it defaults to 2017


testFile=""
# Input source
if options.runOn2017:
    if options.runOnMC:
        testFile = '/store/mc/RunIIFall17MiniAODv2/bbDM_2HDMa_LO_5f_TuneCP3_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_rp_94X_mc2017_realistic_v14-v1/120000/00EC94C6-59F9-EA11-A069-000D3AB6FB6B.root'
    else:
        testFile='/store/data/Run2017B/MET/MINIAOD/31Mar2018-v1/100000/16963797-0937-E811-ABE2-008CFAE45134.root'
elif options.runOn2016:
    if options.runOnMC:
        testFile='/store/mc/RunIIFall17MiniAODv2/bbDM_2HDMa_LO_5f_TuneCP3_13TeV_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_rp_94X_mc2017_realistic_v14-v1/120000/00EC94C6-59F9-EA11-A069-000D3AB6FB6B.root'
    else:
        testFile='/store/data/Run2016B/MET/MINIAOD/17Jul2018_ver2-v1/40000/FE78E8B0-288C-E811-81FC-0025904CDDF8.root'


process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            #fileNames = cms.untracked.vstring("file:/tmp/khurana/temp.root"),
							fileNames = cms.untracked.vstring(testFile),
							#skipEvents = cms.untracked.uint32(0)
                            )


from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJetsAK8'), ## output will be selectedUpdatedPatJets
   pvSource = cms.InputTag('offlineSlimmedPrimaryVertices'),
   svSource = cms.InputTag('slimmedSecondaryVertices'),
   rParam = 0.8,
   #jetCorrections = ('AK8PFchs', cms.vstring(['L2Relative', 'L3Absolute']), 'None'),
   jetCorrections = ('AK8PFPuppi', cms.vstring(['L2Relative', 'L3Absolute','L2L3Residual']), 'None'),
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

      ## for DeepAK8
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight',
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight',
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD',
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD',
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD',
      'pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD',
      ]
        )
## This is for modified MET, needed only for 2017 data
if options.runOn2017:

    process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
    baddetEcallist = cms.vuint32([872439604,872422825,872420274,872423218,
     872423215,872416066,872435036,872439336,
     872420273,872436907,872420147,872439731,
     872436657,872420397,872439732,872439339,
     872439603,872422436,872439861,872437051,
     872437052,872420649,872422436,872421950,
     872437185,872422564,872421566,872421695,
     872421955,872421567,872437184,872421951,
     872421694,872437056,872437057,872437313])

    
    process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
      EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
      ecalMinEt        = cms.double(50.),
      baddetEcal    = baddetEcallist, 
      taggingMode = cms.bool(True),
      debug = cms.bool(False)
      )
    
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
### do we still need this? I guess no.
jetToolbox( process, 'ca15', 'jetSequence', 'out', PUMethod='Puppi', miniAOD=options.useMiniAOD, runOnMC=options.runOnMC,
	    bTagDiscriminators=(bTagDiscriminators + ([] if NOTADDHBBTag else ['pfBoostedDoubleSecondaryVertexCA15BJetTags'])),
	    JETCorrPayload='AK8PFPuppi',JETCorrLevels=jetCorrectionLevelsPuppi,
	    subJETCorrPayload='AK4PFPuppi',subJETCorrLevels=jetCorrectionLevelsPuppi,
	    Cut='pt>120',
	    addSoftDrop=True,addSoftDropSubjets=True, betaCut=1.0, zCutSD=0.15,
	    addNsub=True )


## Jet Energy Resolution
process.patSmearedJets = cms.EDProducer("SmearedPATJetProducer",
    src = cms.InputTag("appliedRegJets"),

    enabled = cms.bool(True),  # If False, no smearing is performed

    rho = cms.InputTag("fixedGridRhoFastjetAll"),

    skipGenMatching = cms.bool(False),  # If True, always skip gen jet matching and smear jet with a random gaussian

    # Resolution and scale factors source.
    # Can be either from GT or text files
    # For GT: only 'algo' must be set
    # For text files: both 'resolutionFile' and 'scaleFactorFile' must point to valid files

    # Read from GT
    algopt = cms.string('AK4PFchs_pt'),
    algo = cms.string('AK4PFchs'),

    # Or from text files
    #resolutionFile = cms.FileInPath('path/to/resolution_file.txt'),
    #scaleFactorFile = cms.FileInPath('path/to/scale_factor_file.txt'),

    # Gen jet matching
    genJets = cms.InputTag("slimmedGenJets"),
    dRMax = cms.double(0.2),  # = cone size (0.4) / 2
    dPtMaxFactor = cms.double(3),  # dPt < 3 * resolution

    # Systematic variation
    # 0: Nominal
    # -1: -1 sigma (down variation)
    # 1: +1 sigma (up variation)
    variation = cms.int32(0),  # If not specified, default to 0

    seed = cms.uint32(37428479),  # If not specified, default to 37428479
    useDeterministicSeed = cms.bool(True),

    debug = cms.untracked.bool(False)
)


## Tau ID embedding

from ExoPieElement.TreeMaker.runTauIdMVA import *
na = TauIDEmbedder(process, cms,
		   debug=True,
		   toKeep = ["2017v2"]
		   )
na.runTauID()

## adding payloads for Tau ID discriminator

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
	src = cms.InputTag("patSmearedJets"),
	levels = jetCorrectionLevelsFullCHS,
	payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJECAK4 = updatedPatJets.clone(
	jetSource = cms.InputTag("patSmearedJets"),
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


process.load('ExoPieElement.TreeMaker.TreeMaker_cfi')
process.tree.useJECText            = cms.bool(options.useJECText)
process.tree.runOn2017             = cms.bool(options.runOn2017)
process.tree.runOn2016             = cms.bool(options.runOn2016)
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
process.tree.fillCA15PuppiJetInfo  = cms.bool(False)

process.tree.THINJets = cms.InputTag("appliedRegJets")
process.tree.FATJets       = cms.InputTag("selectedUpdatedPatJets")#("slimmedJetsAK8")
process.tree.FATJetsForPrunedMass       = cms.InputTag("slimmedJetsAK8")
process.tree.AK4PuppiJets  = cms.InputTag("slimmedJetsPuppi")

if options.runOn2016:
    process.tree.pfType1Met = cms.InputTag("slimmedMETs")
if options.runOn2017:
    if options.runOnMC:
        process.tree.filterLabel = cms.InputTag("TriggerResults::PAT")
    else:
        process.tree.filterLabel = cms.InputTag("TriggerResults::RECO")
    process.tree.pfType1Met = cms.InputTag("slimmedMETsModifiedMET")

# if options.useJECText:
# 	process.tree.THINJets      = cms.InputTag("patSmearedJets")
# 	process.tree.FATJets       = cms.InputTag("selectedUpdatedPatJets")#("slimmedJetsAK8")
# 	process.tree.FATJetsForPrunedMass       = cms.InputTag("slimmedJetsAK8")
# 	process.tree.AK4PuppiJets  = cms.InputTag("slimmedJetsPuppi")




## output file name
process.TFileService = cms.Service("TFileService",fileName = cms.string("ExoPieElementTuples.root"))


##Trigger Filter
if options.runOn2017:
    from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2017BtoF"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)

    process.trigFilter = cms.EDFilter('TrigFilter',
                                      TrigTag = cms.InputTag("TriggerResults::HLT"),
                                      TrigPaths = cms.vstring("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                                                              "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                                              "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight",
                                                              "HLT_Ele27_WPTight_Gsf",
                                                              "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                                                              "HLT_Ele32_WPTight_Gsf",
                                                              "HLT_Ele35_WPTight_Gsf",
                                                              "HLT_IsoMu24",
                                                              "HLT_IsoMu27",
                                                              "HLT_IsoTkMu27",
                                                              "HLT_IsoTkMu24",
                                                              "HLT_Photon200" ),
                                      isMC_ = cms.bool(options.runOnMC)
                                     )
elif options.runOn2016:
    from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
    process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False)
    process.trigFilter = cms.EDFilter('TrigFilter',
                                      TrigTag = cms.InputTag("TriggerResults::HLT"),
                                      TrigPaths = cms.vstring("HLT_PFMET170_BeamHaloCleaned",
                                                              "HLT_PFMET170_HBHE_BeamHaloCleaned",
                                                              "HLT_PFMET170_NotCleaned",
                                                              "HLT_PFMET170_NoiseCleaned",
                                                              "HLT_PFMET170_JetIdCleaned",
                                                              "HLT_PFMET170_HBHECleaned",
                                                              "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight",
                                                              "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight",
                                                              "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight",
                                                              "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",
                                                              "HLT_PFMET110_PFMHT110_IDTight",
                                                              "HLT_IsoMu24","HLT_IsoTkMu24",
                                                              "HLT_IsoMu27","HLT_IsoTkMu27",
                                                              "HLT_Ele27_WPTight_Gsf",
                                                              "HLT_Ele105_CaloIdVT_GsfTrkIdT",
                                                              "HLT_Ele115_CaloIdVT_GsfTrkIdT",
                                                              "HLT_Ele32_WPTight_Gsf",
                                                              "HLT_IsoMu20",
                                                              "HLT_Ele27_eta2p1_WPTight_Gsf",
                                                              "HLT_Ele27_WPLoose_Gsf",
                                                              "HLT_Ele32_eta2p1_WPTight_Gsf",
                                                              "HLT_Photon165_HE10",
                                                              "HLT_Photon175",),
                                      isMC_ = cms.bool(options.runOnMC)
                                     )

process.appliedRegJets= cms.EDProducer('bRegressionProducer',
                                           JetTag=cms.InputTag("slimmedJets"),
                                           rhoFixedGridCollection = cms.InputTag('fixedGridRhoFastjetAll'),
                                           #bRegressionWeightfile= cms.untracked.string("/afs/cern.ch/work/d/dekumar/public/flashgg_setup/CMSSW_8_0_28/src/flashgg/MetaData/data/DNN_models/model-18"),
                                           y_mean = cms.untracked.double(1.0454729795455933) ,
                                           y_std = cms.untracked.double( 0.31628304719924927)
                                           )

if options.runOn2017:
	if not options.useJECText:
		process.analysis = cms.Path(
			process.ecalBadCalibReducedMINIAODFilter*
			process.trigFilter
            *process.prefiringweight
			*process.rerunMvaIsolationSequence
			*process.NewTauIDsEmbedded+
			process.egammaPostRecoSeq+
			process.appliedRegJets+
			process.fullPatMetSequenceModifiedMET+
			#process.patSmearedJets+
			# process.pfMet+
			process.tree
			)
	else:
		process.analysis = cms.Path(
			process.ecalBadCalibReducedMINIAODFilter*
			process.trigFilter
            *process.prefiringweight
			*process.rerunMvaIsolationSequence
			*process.NewTauIDsEmbedded+
			process.egammaPostRecoSeq+
			process.appliedRegJets+
			process.fullPatMetSequenceModifiedMET+
			process.patSmearedJets+
			process.jetCorrSequenceAK4+  ## only when using JEC text files
			process.jetCorrSequenceAK8+  ## only when using JEC text files
			process.jetCorrSequenceAK4Puppi+ ## only when using JEC text files
			process.jetCorrSequenceForPrunedMass+ ## only when using JEC text files
			# process.pfMet+
			process.tree
			)

elif options.runOn2016:
	if not options.useJECText:
		process.analysis = cms.Path(
			process.trigFilter
            *process.prefiringweight
			*process.rerunMvaIsolationSequence
			*process.NewTauIDsEmbedded+
			process.egammaPostRecoSeq+
			process.appliedRegJets+
			#process.patSmearedJets+
			# process.pfMet+
			process.tree
			)
	else:
		process.analysis = cms.Path(
			process.trigFilter
            *process.prefiringweight
			*process.rerunMvaIsolationSequence
			*process.NewTauIDsEmbedded+
			process.egammaPostRecoSeq+
			process.appliedRegJets+
			process.patSmearedJets+
			process.jetCorrSequenceAK4+  ## only when using JEC text files
			process.jetCorrSequenceAK8+  ## only when using JEC text files
			process.jetCorrSequenceAK4Puppi+ ## only when using JEC text files
			process.jetCorrSequenceForPrunedMass+ ## only when using JEC text files
			# process.pfMet+
			process.tree
			)

#print process.dumpPython()
