import FWCore.ParameterSet.Config as cms

process = cms.Process('MVAMET')
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = cms.string("PHYS14_25_V2")


option = 'RECO' # 'GEN' or 'RECO'

### GEN level studies
if option == 'GEN':
    process.load("ExoDiBosonResonances.EDBRGenStudies.genMuons_cff")
    process.load("ExoDiBosonResonances.EDBRGenStudies.genElectrons_cff")
    process.load("ExoDiBosonResonances.EDBRGenStudies.genFatJets_cff")

### RECO level studies
if option == 'RECO':
    process.load("ExoDiBosonResonances.EDBRCommon.goodMuons_cff")
    process.load("ExoDiBosonResonances.EDBRCommon.goodElectrons_cff")
    process.load("ExoDiBosonResonances.EDBRJets.redoSubstructure_cff")
    process.load("ExoDiBosonResonances.EDBRJets.redoPatJets_cff")
    process.load("ExoDiBosonResonances.EDBRCommon.goodJets_cff")
    process.load("ExoDiBosonResonances.EDBRCommon.goodAK4Jets_cff")


# Updates
if option == 'RECO':
    process.goodMuons.src = "slimmedMuons"
    process.goodElectrons.src = "isolatedElectrons"
    process.goodJets.src = "slimmedJetsAK8"
    process.goodAK4Jets.src = "slimmedJets"

process.leptonSequence = cms.Sequence(process.muSequence +
                                      process.eleSequence)



process.jetSequence = cms.Sequence(process.fatJetsSequence + process.substructureSequence + process.redoPatJets + process.ak4JetsSequence)




CA8jecUnc='jec/START53_V23_Uncertainty_AK7PFchs.txt'
CA8jecLevels = [
    'jec/START53_V23_L1FastJet_AK7PFchs.txt',
    'jec/START53_V23_L2Relative_AK7PFchs.txt',
    'jec/START53_V23_L3Absolute_AK7PFchs.txt'
    ]
AK5jecUnc='jec/START53_V23_Uncertainty_AK5PFchs.txt'
AK5jecLevels = [
    'jec/START53_V23_L1FastJet_AK5PFchs.txt',
    'jec/START53_V23_L2Relative_AK5PFchs.txt',
    'jec/START53_V23_L3Absolute_AK5PFchs.txt'
    ]


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)



# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(
        #'file:ZPrimeToTTBarMiniAOD.root'
#        'file:/afs/cern.ch/work/k/khurana/RunIIExoPieElement/CMSSW_7_2_3_patch1/DYJetsMiniAOD/DYJets_100-200_1.root'
        # 'file:0432E62A-7A6C-E411-87BB-002590DB92A8.root' 
        #'file:58462360-8C07-E411-BCE4-008CFA007B98.root'
        #'file:0C3E1051-ED71-E411-B8C9-002590DB91CE.root'
       #'/store/relval/CMSSW_7_4_0/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V7_GENSIM_7_1_15-v1/00000/C0B7E3DE-6FDD-E411-A127-0025905B85EE.root' 
       #'file:C0B7E3DE-6FDD-E411-A127-0025905B85EE.root' 
      # '/store/relval/CMSSW_7_4_0_pre7/RelValADDMonoJet_d3MD3_13/MINIAODSIM/MCRUN2_74_V7-v1/00000/10661B83-43B7-E411-ADEA-0025905B8582.root'
#       '/store/relval/CMSSW_7_4_0/RelValProdTTbar_13/MINIAODSIM/MCRUN2_74_V7_GENSIM_7_1_15-v1/00000/1E12B842-93DD-E411-AF4F-0025905A48D0.root'
       #'/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'
#     '/store/relval/CMSSW_7_4_1/RelValTTbar_13_unsch/MINIAODSIM/MCRUN2_74_V9-v1/00000/EA076017-DCEC-E411-BE66-0025905A60F2.root'  
   #'/store/mc/Phys14DR/DYJetsToEEMuMu_M-1400To2300_13TeV-madgraph/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/12B916AF-CE73-E411-8419-002590747DE2.root' 
  '/store/relval/CMSSW_7_4_0/RelValProdQCD_Pt_3000_3500_13/MINIAODSIM/MCRUN2_74_V7_GENSIM_7_1_15-v1/00000/08D1D655-7BDE-E411-8402-0025905A6060.root'   
  ),
                            skipEvents = cms.untracked.uint32(0)         
                            )


process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)



process.load("RecoJets.JetProducers.ak4PFJets_cfi")
process.ak4PFJets.src = cms.InputTag("packedPFCandidates")

from JetMETCorrections.Configuration.DefaultJEC_cff import ak4PFJetsL1FastL2L3

process.load("RecoMET.METPUSubtraction.mvaPFMET_cff")
#process.pfMVAMEt.srcLeptons = cms.VInputTag("slimmedElectrons")
process.pfMVAMEt.srcPFCandidates = cms.InputTag("packedPFCandidates")
process.pfMVAMEt.srcVertices = cms.InputTag("offlineSlimmedPrimaryVertices")

process.puJetIdForPFMVAMEt.jec =  cms.string('AK4PF')
#process.puJetIdForPFMVAMEt.jets = cms.InputTag("ak4PFJets")
process.puJetIdForPFMVAMEt.vertexes = cms.InputTag("offlineSlimmedPrimaryVertices")
process.puJetIdForPFMVAMEt.rho = cms.InputTag("fixedGridRhoFastjetAll")



# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', 'All')

from RecoMET.METProducers.PFMET_cfi import pfMet
process.pfMet = pfMet.clone(src = "packedPFCandidates")
process.pfMet.calculateSignificance = False # this can't be easily implemented on packed PF candidates at the moment



### for adding jet collection

process.packedGenParticlesForJetsNoNu = cms.EDFilter("CandPtrSelector", src = cms.InputTag("packedGenParticles"), cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"))
## Fat jets (Gen and Reco)
from RecoJets.JetProducers.ak4GenJets_cfi import ak4GenJets
process.genJetsNoNu = ak4GenJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu")
)
from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
process.PFJetsCHS = ak4PFJets.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = getattr(process,"ak4PFJets").src ,
    srcPVs = getattr(process,"ak4PFJets").srcPVs ,
    doAreaFastjet = cms.bool(True),
    jetPtMin = cms.double(150.0)
)
## Pruned fat jets (Gen and Reco) (each module produces two jet collections, fat jets and subjets)
from RecoJets.JetProducers.SubJetParameters_cfi import SubJetParameters
process.genJetsNoNuPruned = ak4GenJets.clone(
    SubJetParameters,
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = cms.InputTag("packedGenParticlesForJetsNoNu") ,
    usePruning = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets")
)
from RecoJets.JetProducers.ak4PFJetsPruned_cfi import ak4PFJetsPruned
process.PFJetsCHSPruned = ak4PFJetsPruned.clone(
    jetAlgorithm = cms.string('AntiKt'),
    rParam = cms.double(0.8),
    src = getattr(process,"ak4PFJets").src ,
    srcPVs = getattr(process,"ak4PFJets").srcPVs ,
    doAreaFastjet = cms.bool(True),
    writeCompound = cms.bool(True),
    jetCollInstanceName=cms.string("SubJets"),
    jetPtMin = cms.double(150.0)
)

postfix = "PFlow"
pfCandidates = 'packedPFCandidates'
pvSource = 'offlineSlimmedPrimaryVertices'
svSource = 'slimmedSecondaryVertices'
muSource = 'slimmedMuons'
elSource = 'slimmedElectrons'
genParticles = 'prunedGenParticles'


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
jetCorrectionsAK4 = ('AK4PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')
jetCorrectionsAK8 = ('AK8PFchs', ['L1FastJet', 'L2Relative', 'L3Absolute'], 'None')


from PhysicsTools.PatAlgos.tools.jetTools import *

addJetCollection(
        process,
        labelName='PFCHS',
        jetSource=cms.InputTag('PFJetsCHS'),
        algo='AK',           # needed for jet flavor clustering
        rParam=0.8, # needed for jet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = False,
        runIVF = False,
        postfix = postfix
    )

getattr(process,'selectedPatJetsPFCHS'+postfix).cut = cms.string('abs(eta) < 2.5')

addJetCollection(
        process,
        labelName='PrunedPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned'),
        algo='AK',
        btagInfos=['None'],
        btagDiscriminators=['None'],
        jetCorrections=jetCorrectionsAK8,
        genJetCollection = cms.InputTag('genJetsNoNu'),
        genParticles = cms.InputTag(genParticles),
        getJetMCFlavour = False, # jet flavor disabled
        postfix = postfix
    )
addJetCollection(
        process,
        labelName='PrunedSubjetsPFCHS',
        jetSource=cms.InputTag('PFJetsCHSPruned','SubJets'),
        algo='AK',           # needed for subjet flavor clustering
        rParam=0.8, # needed for subjet flavor clustering
        pfCandidates = cms.InputTag(pfCandidates),
        pvSource = cms.InputTag(pvSource),
        svSource = cms.InputTag(svSource),
        muSource = cms.InputTag(muSource),
        elSource = cms.InputTag(elSource),
        btagInfos = bTagInfos,
        btagDiscriminators = bTagDiscriminators,
        jetCorrections = jetCorrectionsAK4,
        genJetCollection = cms.InputTag('genJetsNoNuPruned','SubJets'),
        genParticles = cms.InputTag(genParticles),
        explicitJTA = True,  # needed for subjet b tagging
        svClustering = True, # needed for subjet b tagging
        fatJets = cms.InputTag('PFJetsCHS'),              # needed for subjet flavor clustering
        groomedFatJets = cms.InputTag('PFJetsCHSPruned'), # needed for subjet flavor clustering
        runIVF = False,
        postfix = postfix
)






















## Establish references between PATified fat jets and subjets using the BoostedJetMerger

process.selectedPatJetsPrunedPFCHSPacked = cms.EDProducer("BoostedJetMerger",
    jetSrc=cms.InputTag("selectedPatJetsPrunedPFCHS"+postfix),
    subjetSrc=cms.InputTag("selectedPatJetsPrunedSubjetsPFCHS"+postfix)
)

## Pack fat jets with subjets
process.packedPatJetsPFCHS = cms.EDProducer("JetSubstructurePacker",
            jetSrc = cms.InputTag('selectedPatJetsPFCHS'+postfix),
            distMax = cms.double(0.8),
            algoTags = cms.VInputTag(),
            algoLabels = cms.vstring(),
            fixDaughters = cms.bool(False)
)


process.packedPatJetsPFCHS.algoTags.append( cms.InputTag('selectedPatJetsPrunedPFCHSPacked') )
process.packedPatJetsPFCHS.algoLabels.append( 'Pruned' )


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














#process.addjetSequence = cms.Sequence(process.selectedPatJetsPFCHSPacked+process.packedPatJetsPFCHS)




###end of add jet collection

























process.tree = cms.EDAnalyzer(
    'TreeMaker',
    fillPUweightInfo_ = cms.bool(True),
    fillEventInfo_ = cms.bool(True),
    fillGenInfo_   = cms.bool(True),
    fillMuonInfo_  = cms.bool(True),
    fillElecInfo_  = cms.bool(True),
    fillJetInfo_   = cms.bool(True), ## Raman switched it off
    fillMetInfo_   = cms.bool(True),
    fillTrigInfo_  = cms.bool(True),
    fillPhotInfo_  = cms.bool(False),
    fillTauInfo_   = cms.bool(True),
    
    genPartLabel=cms.InputTag("prunedGenParticles"),
    genJetLabel=cms.InputTag("slimmedGenJets"),
    maxNumGenPar  =  cms.uint32(30),
    patMuons=cms.InputTag("slimmedMuons"),
    patElectrons = cms.InputTag("slimmedElectrons"),
    
    #PrunedJets=cms.InputTag("ak8PFJetsCHSPrunedLinks"),
    #PrunedJets=cms.InputTag("slimmedJetsAK8"),
    
    pvSrc  = cms.InputTag('offlineSlimmedPrimaryVertices'),
    tauLabel_ = cms.untracked.InputTag("slimmedTaus"),
    rhoSrc = cms.InputTag('kt6PFJets','rho'),
    ### CA8Jet
    CA8Jets=cms.InputTag("slimmedJetsAK8"),
    CA8jecPayloadNames = cms.vstring( CA8jecLevels ),
    CA8jecUncName = cms.string(CA8jecUnc),    
    ### AK5Jet
    #AK5Jets=cms.InputTag("ak8PFJetsCHS"),
    AK5Jets=cms.InputTag("slimmedJets"),
    AK5jecPayloadNames = cms.vstring( AK5jecLevels ),
    AK5jecUncName = cms.string(AK5jecUnc),    
    
    patMetRaw=cms.InputTag("pfMet"),
    patMet = cms.InputTag("slimmedMETs"),
    #patMet = cms.InputTag("pfMet"),
    #patMetRaw=cms.InputTag("slimmedMETs"),
    
    AK8AddPY= cms.InputTag("packedPatJetsPFCHS"),
    SubJetsPY= cms.InputTag('selectedPatJetsSoftDropPFCHSPacked','SubJets'),
    outFileName=cms.string('outputFileName.root')
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mvamet.root")
                                   )





process.analysis = cms.Path(process.leptonSequence+
                            process.jetSequence+                            
                            #process.addjetSequence+
                            process.pfMVAMEtSequence+
                            process.pfMet+
                            process.tree
                            )


#print process.dumpPython()
