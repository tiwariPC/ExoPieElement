import FWCore.ParameterSet.Config as cms
## removed cleaning from Exo VV package 
## 

process = cms.Process('MVAMET')
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

options.register ('runOn25ns',
		  True ,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "runOn25ns")
options.register ('useMiniAOD',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "useMiniAOD")

options.register ('useJECText',
		  True,
		  VarParsing.multiplicity.singleton,
		  VarParsing.varType.bool,
		  "useJECText")

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


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
if options.runOnMC:
	if options.runOn25ns:
		# 25-ns global tag
		process.GlobalTag = GlobalTag(process.GlobalTag, '80X_mcRun2_asymptotic_2016_miniAODv2_v1', '')  
	else:
		# 50-ns global tag
		process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_startup_v2', '')   
else:## Data no global tag yet
        process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_Prompt_ICHEP16JEC_v0', '')



## for Filters

if options.runOnMC:
	filterlabel="TriggerResults::PAT"
else:
	## for re-reco
	if options.isReReco:
		filterlabel="TriggerResults::RECO"
	## for prompt-reco
	else:
		filterlabel="TriggerResults::PAT"


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
	testFile='/store/mc/RunIISpring16MiniAODv2/ZprimeToA0hToA0chichihbb_2HDM_MZp-1000_MA0-300_13TeV-madgraph/MINIAODSIM/PUSpring16RAWAODSIM_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/10000/4E2DAE88-FE25-E611-BABE-02163E01186C.root'
else:
	testFile='/store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v2/000/273/158/00000/02D9C19F-571A-E611-AD8E-02163E013732.root'

process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(testFile),
			    #skipEvents = cms.untracked.uint32(0)         
                            )




## skip the events 
## using the MET tails events
if options.runOnMC:
    print "No events to skip"
else:
    rangeEventsToSkip = cms.untracked.VEventRange(listEventsToSkip)
    process.source.eventsToSkip = rangeEventsToSkip

## For MVA MET
## Check this reciepe when using MVA MET
process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.ak4PFJets.src = cms.InputTag("packedPFCandidates")

from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
#process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
process.pfMVAMEt.inputFileNames = cms.PSet(
        U     = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru_7_4_X_miniAOD_25NS_July2015.root'),
        DPhi  = cms.FileInPath('RecoMET/METPUSubtraction/data/gbrphi_7_4_X_miniAOD_25NS_July2015.root'),
        CovU1 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru1cov_7_4_X_miniAOD_25NS_July2015.root'),
        CovU2 = cms.FileInPath('RecoMET/METPUSubtraction/data/gbru2cov_7_4_X_miniAOD_25NS_July2015.root')
	)

process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
#process.puJetIdForPFMVAMEt.jets = cms.InputTag("ak4PFJets")
process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")


## MVA MET Ends Here 
##



##
## This is for Uncorrected MET
from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment
## Uncorrected MET edns here 
##




'''

### for adding jet collection

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## produce Fat ak8 jets (Gen and Reco)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.genJetsNoNuAK8 = ak4GenJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.PFJetsCHSAK8 = ak4PFJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = getattr(process,"ak4PFJets").src ,
    srcPVs = getattr(process,"ak4PFJets").srcPVs ,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(150.0)
)
## produce Pruned ak8 fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPrunedAK8 = ak4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu") ,
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.PFJetsCHSPrunedAK8 = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = getattr(process,"ak4PFJets").src ,
    srcPVs = getattr(process,"ak4PFJets").srcPVs ,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(150.0)
)


## Produce SoftDrop ak8 fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
process.genJetsNoNuSoftDrop = ak4GenJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu"),
    useSoftDrop = cms.bool(False),# because it crashes for PHYS14 samples
    zcut = cms.double(0.1),
    beta = cms.double(0.0),
    R0 = cms.double(0.8),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsSoftDrop_cfi import ak4PFJetsSoftDrop
process.PFJetsCHSSoftDrop = ak4PFJetsSoftDrop.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    R0 = cms.double(0.8),
    src = getattr(process,"ak4PFJets").src,
    srcPVs = getattr(process,"ak4PFJets").srcPVs,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(150.0)
)











'''

postfix = "PFlow"
pfCandidates = 'packedPFCandidates'
jetSource = 'ak4PFJets'
pvSource = 'offlineSlimmedPrimaryVertices'
svSource = 'slimmedSecondaryVertices'
muSource = 'slimmedMuons'
elSource = 'slimmedElectrons'
genParticles = 'prunedGenParticles'
genJetCollection = 'ak4GenJetsNoNu'


bTagInfos = [
    'pfImpactParameterTagInfos'
   ,'pfSecondaryVertexTagInfos'
   ,'pfInclusiveSecondaryVertexFinderTagInfos'
   ,'softPFMuonsTagInfos'
   ,'softPFElectronsTagInfos'
]
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
]

## Jet energy corrections

## For jet energy correction
if options.runOnMC:
	jetCorrectionsAK4CHS       = ('AK4PFchs', ['L1FastJet','L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8CHS       = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8CHSL23    = ('AK8PFchs', ['L2Relative', 'L3Absolute'], 'None')
	jetCorrectionsAK8PuppiL23  = ('AK8PFPuppi', ['L2Relative', 'L3Absolute'], 'None')
	jetCorrectionLevelsFullCHS = ['L1FastJet', 'L2Relative', 'L3Absolute']
	jetCorrectionLevels23CHS   = ['L2Relative', 'L3Absolute']

	AK4JECTextFiles = [
		'Spring16_25nsV6_MC_L1FastJet_AK4PFchs.txt',
		'Spring16_25nsV6_MC_L2Relative_AK4PFchs.txt',
		'Spring16_25nsV6_MC_L3Absolute_AK4PFchs.txt'
		]
	AK4JECUncTextFile = 'Spring16_25nsV6_MC_Uncertainty_AK4PFchs.txt'

	AK8JECTextFiles = [
		'Spring16_25nsV6_MC_L1FastJet_AK8PFchs.txt',
		'Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt',
		'Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt'
		]
	AK8JECUncTextFile = 'Spring16_25nsV6_MC_Uncertainty_AK8PFchs.txt'  
	prunedMassJECTextFiles = [
		'Spring16_25nsV6_MC_L2Relative_AK8PFchs.txt',
		'Spring16_25nsV6_MC_L3Absolute_AK8PFchs.txt'
		]
	softdropMassJECTextFiles = [
		'Spring16_25nsV6_MC_L2Relative_AK8PFPuppi.txt',
		'Spring16_25nsV6_MC_L3Absolute_AK8PFPuppi.txt'
		]

else:
        jetCorrectionsAK4CHS       = ('AK4PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8CHS       = ('AK8PFchs', ['L1FastJet','L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8CHSL23    = ('AK8PFchs', ['L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionsAK8PuppiL23  = ('AK8PFPuppi', ['L2Relative', 'L3Absolute','L2L3Residual'], 'None')
	jetCorrectionLevelsFullCHS = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']
	jetCorrectionLevels23CHS   = ['L2Relative', 'L3Absolute','L2L3Residual']
	AK4JECTextFiles = [
		'Spring16_25nsV6_DATA_L1FastJet_AK4PFchs.txt',
		'Spring16_25nsV6_DATA_L2Relative_AK4PFchs.txt',
		'Spring16_25nsV6_DATA_L3Absolute_AK4PFchs.txt',
		'Spring16_25nsV6_DATA_L2L3Residual_AK4PFchs.txt'
		]
	AK4JECUncTextFile = 'Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt'
	AK8JECTextFiles = [
		'Spring16_25nsV6_DATA_L1FastJet_AK8PFchs.txt',
		'Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt',
		'Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt',
		'Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt'
		]
	AK8JECUncTextFile = 'Spring16_25nsV6_DATA_Uncertainty_AK8PFchs.txt'
	prunedMassJECTextFiles = [
		'Spring16_25nsV6_DATA_L2Relative_AK8PFchs.txt',
		'Spring16_25nsV6_DATA_L3Absolute_AK8PFchs.txt',
		'Spring16_25nsV6_DATA_L2L3Residual_AK8PFchs.txt'
		]
	softdropMassJECTextFiles = [
		'Spring16_25nsV6_DATA_L2Relative_AK8PFPuppi.txt',
		'Spring16_25nsV6_DATA_L3Absolute_AK8PFPuppi.txt',
		'Spring16_25nsV6_DATA_L2L3Residual_AK8PFPuppi.txt'
		]


from PhysicsTools.PatAlgos.tools.jetTools import *


NOTADDHBBTag=False
'''
# add PFJetsCHSAK8

addJetCollection(
        process,
        labelName='PFCHSAK8',
        jetSource=cms.InputTag('PFJetsCHSAK8'),
        algo='AK',           # needed for jet flavor clustering
        rParam=0.8, # needed for jet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        #btagDiscriminators = bTagDiscriminators,
        btagDiscriminators = (bTagDiscriminators + ([] if NOTADDHBBTag else ['pfBoostedDoubleSecondaryVertexAK8BJetTags'])),
        jetCorrections = jetCorrectionsAK8CHS,
        genJetCollection = cms.InputTag('genJetsNoNuAK8'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = False,
        runIVF = False,
        postfix = postfix
    )

# used pat selector on PFCHSAK8
getattr(process,'selectedPatJetsPFCHSAK8'+postfix).cut = cms.string('abs(eta) < 2.5')


# Add SoftDrop jet
addJetCollection(
        process,
        labelName='SoftDropPFCHS',
        jetSource=cms.InputTag('PFJetsCHSSoftDrop'),
        algo='AK',
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8CHS,
        genJetCollection = cms.InputTag('genJetsNoNuAK8'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )

# Add SoftDrop subjet
addJetCollection(
        process,
        labelName='SoftDropSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSSoftDrop','SubJets'),
        algo='AK',           # needed for subjet flavor clustering
        rParam=0.8, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4CHS,
        genJetCollection = cms.InputTag('genJetsNoNuSoftDrop','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHSAK8'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSSoftDrop'), # needed for subjet flavor clustering
        runIVF = False,
        postfix = postfix
    )

## Softdrop  Establish references between PATified fat jets and subjets using the BoostedJetMerger
process.selectedPatJetsSoftDropPFCHSPacked = cms.EDProducer("BoostedJetMerger",
        jetSrc=cms.InputTag("selectedPatJetsSoftDropPFCHS"+postfix),
        subjetSrc=cms.InputTag("selectedPatJetsSoftDropSubjetsPFCHS"+postfix)
)



#add pruned PFJetsCHSAK8
addJetCollection(
        process,
        labelName='PrunedPFCHSAK8',
        jetSource=cms.InputTag('PFJetsCHSPrunedAK8'),
        algo='AK',
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8CHS,
        genJetCollection = cms.InputTag('genJetsNoNuAK8'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )

#add subjet from pruned PFJetsCHSAK8

addJetCollection(
        process,
        labelName='PrunedSubjetsPFCHSAK8',
        jetSource=cms.InputTag('PFJetsCHSPrunedAK8','SubJets'),
        algo='AK',           # needed for subjet flavor clustering
        rParam=0.8, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4CHS,
        genJetCollection = cms.InputTag('genJetsNoNuPrunedAK8','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHSAK8'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSPrunedAK8'), # needed for subjet flavor clustering
        runIVF = False,
        postfix = postfix
)






## Establish references between PATified fat jets and subjets using the BoostedJetMerger

process.selectedPatJetsPrunedPFCHSAK8Packed = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHSAK8"+postfix),
    subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHSAK8"+postfix)
)



## Pack fat jets with subjets
process.packedPatJetsPFCHSAK8 = cms.EDProducer("JetSubstructurePacker",
            jetSrc = cms.InputTag('selectedPatJetsPFCHSAK8'+postfix),
            distMax = cms.double(0.8),
            algoTags = cms.VInputTag(),
            algoLabels = cms.vstring(),
            fixDaughters = cms.bool(False)
)

SaveSoftDrop=True;

if SaveSoftDrop:
        process.packedPatJetsPFCHSAK8.algoTags.append( cms.InputTag('selectedPatJetsSoftDropPFCHSPacked') )
        process.packedPatJetsPFCHSAK8.algoLabels.append( 'SoftDrop' )
else:
        process.packedPatJetsPFCHSAK8.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSAK8Packed') )
        process.packedPatJetsPFCHSAK8.algoLabels.append( 'Pruned' )



#process.packedPatJetsPFCHSAK8.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSAK8Packed') )
#process.packedPatJetsPFCHSAK8.algoLabels.append( 'Pruned' )


for m in ['patJets'+postfix, 'patJetsPFCHSAK8'+postfix, 'patJetsSoftDropSubjetsPFCHS'+postfix, 'patJetsPrunedSubjetsPFCHSAK8'+postfix]:
    if hasattr(process,m) and getattr( getattr(process,m), 'addBTagInfo' ):
        print "Switching 'addTagInfos' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addTagInfos', cms.bool(True) )
    if hasattr(process,m):
        print "Switching 'addJetFlavourInfo' for " + m + " to 'True'"
        setattr( getattr(process,m), 'addJetFlavourInfo', cms.bool(True) )



jetRadius=0.8

        # Set the cone size for the jet-track association to the jet radius
getattr(process,'pfImpactParameterTagInfosPFCHSAK8'+postfix).maxDeltaR = cms.double(jetRadius) # default is 0.4
getattr(process,'pfSecondaryVertexTagInfosPFCHSAK8'+postfix).trackSelection.jetDeltaRMax = cms.double(jetRadius)   # default is 0.3
getattr(process,'pfSecondaryVertexTagInfosPFCHSAK8'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(jetRadius) # default is 0.4
        # Set the jet-SV dR to the jet radius
getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHSAK8'+postfix).vertexCuts.maxDeltaRToJetAxis = cms.double(jetRadius) # default is 0.4
getattr(process,'pfInclusiveSecondaryVertexFinderTagInfosPFCHSAK8'+postfix).extSVDeltaRToJet = cms.double(jetRadius) # default is 0.3
        # Set the JP track dR cut to the jet radius
process.candidateJetProbabilityComputerFat = process.candidateJetProbabilityComputer.clone( deltaR = cms.double(jetRadius) ) # default is 0.3
getattr(process,'pfJetProbabilityBJetTagsPFCHSAK8'+postfix).jetTagComputer = cms.string('candidateJetProbabilityComputerFat')
        # Set the JBP track dR cut to the jet radius
process.candidateJetBProbabilityComputerFat = process.candidateJetBProbabilityComputer.clone( deltaR = cms.double(jetRadius) ) # default is 0.4
getattr(process,'pfJetBProbabilityBJetTagsPFCHSAK8'+postfix).jetTagComputer = cms.string('candidateJetBProbabilityComputerFat')
        # Set the CSV track dR cut to the jet radius
process.candidateCombinedSecondaryVertexComputerFat = process.candidateCombinedSecondaryVertexComputer.clone()
process.candidateCombinedSecondaryVertexComputerFat.trackSelection.jetDeltaRMax = cms.double(jetRadius) # default is 0.3
process.candidateCombinedSecondaryVertexComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(jetRadius) # default is 0.3
getattr(process,'pfCombinedSecondaryVertexV2BJetTagsPFCHSAK8'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexComputerFat')
        # Set the CSVv2 track dR cut to the jet radius
process.candidateCombinedSecondaryVertexV2ComputerFat = process.candidateCombinedSecondaryVertexV2Computer.clone()
process.candidateCombinedSecondaryVertexV2ComputerFat.trackSelection.jetDeltaRMax = cms.double(jetRadius) # default is 0.3
process.candidateCombinedSecondaryVertexV2ComputerFat.trackPseudoSelection.jetDeltaRMax = cms.double(jetRadius) # default is 0.3
getattr(process,'pfCombinedInclusiveSecondaryVertexV2BJetTagsPFCHSAK8'+postfix).jetTagComputer = cms.string('candidateCombinedSecondaryVertexV2ComputerFat')




from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

process.goodaddJets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                        filterParams = pfJetIDSelector.clone(),
                        src = cms.InputTag("packedPatJetsPFCHSAK8")
                        )

### Cleaning
# We want to make sure that the jets are not the electrons or muons done previously

import PhysicsTools.PatAlgos.cleaningLayer1.jetCleaner_cfi as jetCleaner_cfi

process.cleanaddJets = jetCleaner_cfi.cleanPatJets.clone()
process.cleanaddJets.src = "goodaddJets"
process.cleanaddJets.checkOverlaps.muons.src = "slimmedMuons"
process.cleanaddJets.checkOverlaps.muons.deltaR = 0.8
process.cleanaddJets.checkOverlaps.muons.requireNoOverlaps = True
process.cleanaddJets.checkOverlaps.electrons.src = "slimmedElectrons"
process.cleanaddJets.checkOverlaps.electrons.deltaR = 0.8
process.cleanaddJets.checkOverlaps.electrons.requireNoOverlaps = True
process.cleanaddJets.checkOverlaps.photons = cms.PSet()
process.cleanaddJets.checkOverlaps.taus = cms.PSet()
process.cleanaddJets.checkOverlaps.tkIsoElectrons = cms.PSet()
process.cleanaddJets.finalCut = "pt > 20 & abs(eta) < 2.4"

process.fataddJetsSequence = cms.Sequence(
               process.goodaddJets + 
               process.cleanaddJets
               )

















#process.patJetsAK8.addTagInfos = cms.bool(True)


'''

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


'''
addJetCollection(
        process,
        labelName='PFCHSAK8',
        jetSource=cms.InputTag('PFJetsCHSAK8'),
        algo='AK',           # needed for jet flavor clustering
        rParam=0.8, # needed for jet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        #btagDiscriminators = bTagDiscriminators,
        btagDiscriminators = (bTagDiscriminators + ([] if NOTADDHBBTag else ['pfBoostedDoubleSecondaryVertexAK8BJetTags'])),
        jetCorrections = jetCorrectionsAK8CHS,
        genJetCollection = cms.InputTag('genJetsNoNuAK8'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = False,
        runIVF = False,
        postfix = postfix
    )


process.miniAODjetSequence = cms.Sequence(
                             process.selectedPatJetsPFCHSAK8PFlow+
                             process.selectedPatJetsPrunedPFCHSAK8Packed +
                             process.fataddJetsSequence
                             )




###end of add jet collection
### FOR adding new hbb b-tags in miniaod?

#process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
#process.load("PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff")

#from RecoJets.Configuration.RecoPFJets_cff import ak8PFJetsCHS



'''

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
if options.runOn25ns:
	my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
			 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
			 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_nonTrig_V1_cff']
else:
	my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_50ns_V2_cff',
			 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff',
			 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_50ns_Trig_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

if options.runOn25ns:
    #my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_PHYS14_PU20bx25_V2_cff']
    my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff',
                        'RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']
else:
    my_phoid_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_50ns_V1_cff']

#add them to the VID producer
for idmod in my_phoid_modules:
	setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)




#process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("ncuslimmedPhoton")
#process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("ncuslimmedElectron")

## For normal AK4 jets jet energy correction on top of miniAOD
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
process.patJetCorrFactorsReapplyJECAK4 = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJets"),
	levels = jetCorrectionLevelsFullCHS,
	payload = 'AK4PFchs' ) # Make sure to choose the appropriate levels and payload here!

from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process.patJetsReapplyJECAK4 = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJets"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK4"))
  )

process.jetCorrSequenceAK4 = cms.Sequence( process.patJetCorrFactorsReapplyJECAK4 + process.patJetsReapplyJECAK4 )




### For normal AK8 jet energy correction on top of miniAOD
process.patJetCorrFactorsReapplyJECAK8 = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJetsAK8"),
	levels = jetCorrectionLevelsFullCHS,
	payload = 'AK8PFchs' ) # Make sure to choose the appropriate levels and payload here!

process.patJetsReapplyJECAK8 = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJetsAK8"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECAK8"))
  )


process.jetCorrSequenceAK8 = cms.Sequence( process.patJetCorrFactorsReapplyJECAK8 + process.patJetsReapplyJECAK8 )



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


## For correcting softdrop jet mass + Puppi
process.patJetCorrFactorsReapplyJECForSoftDropMass = updatedPatJetCorrFactors.clone(
	src = cms.InputTag("slimmedJetsAK8"),
	levels = jetCorrectionLevels23CHS,
	payload = 'AK8PFPuppi' ) # Make sure to choose the appropriate levels and payload here!


process.patJetsReapplyJECForSoftDropMass = updatedPatJets.clone(
	jetSource = cms.InputTag("slimmedJetsAK8"),
	jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECForSoftDropMass"))
	)

process.jetCorrSequenceForSoftDropMass = cms.Sequence( process.patJetCorrFactorsReapplyJECForSoftDropMass + process.patJetsReapplyJECForSoftDropMass )



process.load('ExoPieElement.TreeMaker.TreeMaker_cfi')
process.tree.filterLabel           = cms.InputTag(filterlabel)
process.tree.useJECText            = cms.bool(options.useJECText)
process.tree.THINjecNames          = cms.vstring(AK4JECTextFiles)
process.tree.THINjecUncName        = cms.string(AK4JECUncTextFile)
process.tree.FATprunedMassJecNames = cms.vstring(prunedMassJECTextFiles)
process.tree.FATsoftdropMassJecNames = cms.vstring(softdropMassJECTextFiles)
process.tree.FATjecNames           = cms.vstring(AK8JECTextFiles)
process.tree.FATjecUncName         = cms.string(AK8JECUncTextFile)
process.tree.ADDjecNames           = cms.vstring(AK8JECTextFiles)
process.tree.ADDjecUncName         = cms.string(AK8JECUncTextFile)

if options.useJECText:
	process.tree.THINJets      = cms.InputTag("slimmedJets")
	process.tree.FATJets       = cms.InputTag("slimmedJetsAK8")
	process.tree.FATJetsForPrunedMass       = cms.InputTag("slimmedJetsAK8")
	process.tree.FATJetsForSoftDropMass     = cms.InputTag("slimmedJetsAK8")

if not options.runOn25ns:
### Electron
	process.tree.eleVetoIdMap      = cms.InputTag  ("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-veto")
	process.tree.eleLooseIdMap     = cms.InputTag ("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-loose")
        process.tree.eleMediumIdMap    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-medium")
        process.tree.eleTightIdMap     = cms.InputTag ("egmGsfElectronIDs:cutBasedElectronID-Spring15-50ns-V2-standalone-tight")
        process.tree.eleMVAMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-50ns-Trig-V1-wp90")
        process.tree.eleMVATightIdMap  = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-50ns-Trig-V1-wp80")
        process.tree.mvaValuesMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig50nsV1Values")
        process.tree.mvaCategoriesMap  = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories")    
### Photon
        process.tree.phoLooseIdMap     = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-loose")
        process.tree.phoMediumIdMap    = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-medium")
        process.tree.phoTightIdMap     = cms.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring15-50ns-V1-standalone-tight")


process.TFileService = cms.Service("TFileService",
				   fileName = cms.string("NCUGlobalTuples.root")          
				   # fileName = cms.string('$outputFileName')          
				   )



process.allEventsCounter = cms.EDFilter(
	"EventCounter"
 )


if not options.useJECText:
	process.analysis = cms.Path(
		process.allEventsCounter+
		#process.ncuslimmer+
		process.egmGsfElectronIDSequence+## by raman
		process.egmPhotonIDSequence+ ## by raman
		#    process.pfMVAMEtSequence+   # disabled before the official code is fixed
		process.pfMet+
		# process.miniAODjetSequence+   ## by raman        
		process.jetCorrSequenceAK4+
		process.jetCorrSequenceAK8+
		process.jetCorrSequenceForPrunedMass+
		process.jetCorrSequenceForSoftDropMass+
		#process.HBHENoiseFilterResultProducer+ ## by raman
		process.tree
		)
else:
	process.analysis = cms.Path(
		process.allEventsCounter+
		process.egmGsfElectronIDSequence+## by raman
		process.egmPhotonIDSequence+ ## by raman
		#    process.pfMVAMEtSequence+   # disabled before the official code is fixed
		process.pfMet+
		# process.miniAODjetSequence+   ## by raman        
		#process.HBHENoiseFilterResultProducer+ ## by raman
		process.tree
		)


#print process.dumpPython()
