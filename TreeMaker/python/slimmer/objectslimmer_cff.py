import FWCore.ParameterSet.Config as cms


import PhysicsTools.PatAlgos.selectionLayer1.jetSelector_cfi as jetSelector_cfi
ncuslimmedAK4Jet = jetSelector_cfi.selectedPatJets.clone()
ncuslimmedAK4Jet.src = cms.InputTag("slimmedJets")
ncuslimmedAK4Jet.cut = cms.string("pt > 20 & abs(eta) < 2.5")


import PhysicsTools.PatAlgos.selectionLayer1.photonSelector_cfi as photonSelector_cfi
ncuslimmedPhoton = photonSelector_cfi.selectedPatPhotons.clone()
ncuslimmedPhoton.src = cms.InputTag("slimmedPhotons")
ncuslimmedPhoton.cut = cms.string("pt > 10. & abs(eta) < 2.7")

import PhysicsTools.PatAlgos.selectionLayer1.electronSelector_cfi as electronSelector_cfi
ncuslimmedElectron = electronSelector_cfi.selectedPatElectrons.clone()
ncuslimmedElectron.src = cms.InputTag("slimmedElectrons")
ncuslimmedElectron.cut = cms.string("pt > 10. & abs(eta) < 2.7")


import PhysicsTools.PatAlgos.selectionLayer1.muonSelector_cfi as muonSelector_cfi
ncuslimmedMuon = muonSelector_cfi.selectedPatMuons.clone()
ncuslimmedMuon.src = cms.InputTag("slimmedMuons")
ncuslimmedMuon.cut = cms.string("pt > 10. & abs(eta) < 2.5")


import PhysicsTools.PatAlgos.selectionLayer1.tauSelector_cfi as tauSelector_cfi
ncuslimmedTau     = tauSelector_cfi.selectedPatTaus.clone()
ncuslimmedTau.src = cms.InputTag("NewTauIDsEmbedded")
ncuslimmedTau.cut = cms.string("pt > 10. & abs(eta) < 2.5")

ncuslimmer = cms.Sequence(ncuslimmedAK4Jet+
                          ncuslimmedPhoton+
                          ncuslimmedElectron+
                          ncuslimmedMuon+
                          ncuslimmedTau)
