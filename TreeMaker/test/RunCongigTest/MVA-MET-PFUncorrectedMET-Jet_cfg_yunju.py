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
    input = cms.untracked.int32(10)
)



# Input source
process.source = cms.Source("PoolSource",
                            secondaryFileNames = cms.untracked.vstring(),
                            fileNames = cms.untracked.vstring(
        #'/store/relval/CMSSW_7_4_0/RelValProdQCD_Pt_3000_3500_13/MINIAODSIM/MCRUN2_74_V7_GENSIM_7_1_15-v1/00000/08D1D655-7BDE-E411-8402-0025905A6060.root'   
        #'file:step3_miniAOD_M1500_9.root'
        'root://eoscms.cern.ch//eos/cms/store/user/khurana/MonoHStep3/step3_miniAOD_M1500_9.root'
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

process.load("ExoPieElement.TreeMaker.TreeMaker_cfi")


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("mvamet.root")
                                   )





process.analysis = cms.Path(process.leptonSequence+
                            process.jetSequence+                            
                            process.pfMVAMEtSequence+
                            process.pfMet+
                            process.tree
                            )


#print process.dumpPython()
