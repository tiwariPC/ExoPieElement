import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
process = cms.Process("OWNJets")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       # 'root://xrootd.cmsaf.mit.edu//store/user/bmaier/moriond17/2HDMa_gg_tb_1p0_MH3_600_MH4_100_MH2_600_MHC_600/2HDMa_gg_tb_1p0_MH3_600_MH4_100_MH2_600_MHC_600_42732727_miniaod.root'
       'root://cms-xrd-global.cern.ch///store/mc/RunIISummer16MiniAODv2/BBbarDMJets_pseudo_NLO_Mchi-1_Mphi-50_TuneCUETP8M1_13TeV-madgraph-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/100000/0470A659-F1CF-E611-BA86-002590E7DFD6.root'
    )
)

#JetCollectionVInputTag      = cms.VInputTag(cms.InputTag('slimmedJets'))


process.myLabel = cms.EDAnalyzer('ProducerTest',
            JetTag      = cms.InputTag('slimmedJets'),
)

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.root')
)


process.p = cms.Path(process.myProducerLabel)

process.e = cms.EndPath(process.out)
