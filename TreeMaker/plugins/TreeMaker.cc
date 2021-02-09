// -*- C++ -*-
//
// Package:    TreeMaker
// Class:      TreeMaker
// Original Author:  Shin-Shan Yu, Raman Khurana, Yun-Ju Lu
//                   National Central University
//         Created:  Tue Jul  6 21:04:59 CEST 2010
//Updated on: July 31 2019, Raman KHURANA


// system include files
#include <memory>
#include <string>
#include <vector>
#include <sstream>

#include "ExoPieElement/TreeMaker/interface/TreeMaker.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
TreeMaker::TreeMaker(const edm::ParameterSet& iConfig):
  genLumiHeaderToken_(consumes<GenLumiInfoHeader,edm::InLumi>(edm::InputTag("generator"))),
  runOnrp2HDM_(iConfig.getParameter<bool>("runOnrp2HDM")),
  runOnrpZpB_(iConfig.getParameter<bool>("runOnrpZpB"))
{
  fillPUweightInfo_=false;
  fillEventInfo_   =false;
  fillMetInfo_     =false;
  fillTrigInfo_    =false;
  fillFilterInfo_  =false;

  fillGenInfo_     =false;
  fillElecInfo_    =false;
  fillMuonInfo_    =false;
  fillTauInfo_     =false;
  fillPhotInfo_    =false;
  fillJetInfo_     =false;
  fillFATJetInfo_  =false;
  fillAK4PuppiJetInfo_ = false;
  fillAK8PuppiJetInfo_ = false;
  fillCA15PuppiJetInfo_ = false;

  fillPUweightInfo_ = iConfig.getParameter<bool>("fillPUweightInfo");
  fillEventInfo_    = iConfig.getParameter<bool>("fillEventInfo");
  fillMetInfo_      = iConfig.getParameter<bool>("fillMetInfo");
  fillTrigInfo_     = iConfig.getParameter<bool>("fillTrigInfo");
  fillFilterInfo_   = iConfig.getParameter<bool>("fillFilterInfo");
  fillGenInfo_      = iConfig.getParameter<bool>("fillGenInfo");
  fillElecInfo_     = iConfig.getParameter<bool>("fillElecInfo");
  fillMuonInfo_     = iConfig.getParameter<bool>("fillMuonInfo");
  fillTauInfo_      = iConfig.getParameter<bool>("fillTauInfo");
  fillPhotInfo_     = iConfig.getParameter<bool>("fillPhotInfo");
  fillJetInfo_      = iConfig.getParameter<bool>("fillJetInfo");
  fillFATJetInfo_   = iConfig.getParameter<bool>("fillFATJetInfo");
  fillAK4PuppiJetInfo_   = iConfig.getParameter<bool>("fillAK4PuppiJetInfo");
  fillAK8PuppiJetInfo_   = iConfig.getParameter<bool>("fillAK8PuppiJetInfo");
  fillCA15PuppiJetInfo_   = iConfig.getParameter<bool>("fillCA15PuppiJetInfo");
  std::cout<<" called ca15 from Treemaker "<<std::endl;

  edm::Service<TFileService> fs;


  bool debug__ = false;
  tree_ = fs->make<TTree>("treeMaker","tree");
  if( fillPUweightInfo_)
    {
      if (debug__) std::cout<< " fillPUweightInfo_ "<<std::endl;
      puweight_                   = new puweight("pu_",tree_);
      puweight_->puInfoToken      = consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"));
    }

  if( fillEventInfo_ )
    {
      if (debug__) std::cout<< " fillEventInfo_"<<std::endl;
      eventInfo_                  = new eventInfo("",tree_);
      eventInfo_->vertexToken     = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
    }


  if( fillMetInfo_ )
    {
      if (debug__) std::cout<< " fillMetInfo_"<<std::endl;
      patMetTree_                 = new patMetTree("pf",tree_);
      patMetTree_->pfMETToken     = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("patMet"));
      patMetTree_->pfMETModifiedToken  = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("pfType1Met"));
      patMetTree_->puppimetToken  = consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("puppiMET"));

    }


  if( fillTrigInfo_ )
    {
      if (debug__) std::cout<< "fillTrigInfo_ "<<std::endl;
      patHltTree_                             = new patHltTree("hlt_",tree_,iConfig);
      patHltTree_->trigResultsToken           = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerLabel"));
      patHltTree_->triggerPrescalesToken      = consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"));
      patHltTree_->triggerObjectsToken         = consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("slimmedPatTrigger"));

    }
  if( fillFilterInfo_ )
    {
      if (debug__) std::cout<< " fillFilterInfo_"<<std::endl;
      patFilterTree_                          = new patFilters("hlt_",tree_);
      patFilterTree_->filterTrigResultsToken  = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("filterLabel"));
      patFilterTree_->ecalBadCalibFilterUpdate_token= consumes< bool >(edm::InputTag("ecalBadCalibReducedMINIAODFilter"));
    }

  if( fillGenInfo_ )
    {
      if (debug__) std::cout<< " fillGenInfo_ "<<std::endl;
      genInfoTree_                           = new genInfoTree("",tree_,iConfig);
      genInfoTree_->genParticleToken         = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genPartLabel"));
      genInfoTree_->genEventToken            = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
      genInfoTree_->lheRunToken              = consumes<LHERunInfoProduct,edm::InRun>(edm::InputTag("externalLHEProducer"));
      genInfoTree_->lheEventToken            = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
      genInfoTree_->genMETToken_true         = consumes<pat::PackedGenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
      genInfoTree_->genMETToken_calo         = consumes<reco::GenMETCollection>(edm::InputTag("genMetCalo"));
      genInfoTree_->genMETToken_caloNonPrompt = consumes<reco::GenMETCollection>(edm::InputTag("genMetCaloAndNonPrompt"));
      genInfoTree_->ak4genJetsToken           = consumes<reco::GenJetCollection>(edm::InputTag("ak4GenJets"));
      genInfoTree_->ak8genJetsToken           = consumes<reco::GenJetCollection>(edm::InputTag("ak8GenJets"));
  }

  if( fillElecInfo_ )
    {
      if (debug__) std::cout<< " fillElecInfo_"<<std::endl;

      patElecTree_                              = new patElecTree("",tree_,iConfig);
      patElecTree_->vertexToken                 = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      patElecTree_->rhoForLepToken              = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
      patElecTree_->eleToken                    = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("eleLabel"));
      patElecTree_->pfCandToken                 = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfForMiniIso"));

      /*
      patElecTree_->eleVetoIdMapToken           = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"));
      patElecTree_->eleLooseIdMapToken          = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"));
      patElecTree_->eleMediumIdMapToken         = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"));
      patElecTree_->eleTightIdMapToken          = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleTightIdMap"));
      patElecTree_->eleHEEPIdMapToken           = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"));

      patElecTree_->eleVetoIdCFToken            = consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("eleVetoIdMap"));
      patElecTree_->eleLooseIdCFToken           = consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("eleLooseIdMap"));
      patElecTree_->eleMediumIdCFToken          = consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("eleMediumIdMap"));
      patElecTree_->eleTightIdCFToken           = consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("eleTightIdMap"));
      patElecTree_->eleHEEPIdCFToken            = consumes<edm::ValueMap<vid::CutFlowResult>>(iConfig.getParameter<edm::InputTag>("eleHEEPIdMap"));

      patElecTree_->eleMVAMediumIdMapToken      = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleMVAMediumIdMap"));
      patElecTree_->eleMVATightIdMapToken       = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("eleMVATightIdMap"));
      patElecTree_->mvaValuesMapToken           = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mvaValuesMap"));
      patElecTree_->mvaCategoriesMapToken       = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("mvaCategoriesMap"));
      */
    }

  if( fillMuonInfo_ )
    {
      if (debug__) std::cout<< " fillMuonInfo_"<<std::endl;
      patMuTree_                             = new patMuonTree("",tree_,iConfig);
      patMuTree_->vertexToken                = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      patMuTree_->rhoForLepToken             = consumes<double>(edm::InputTag("fixedGridRhoFastjetCentralNeutral"));
      patMuTree_->muToken                    = consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muoLabel"));
      patMuTree_->pfCandToken                = consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfForMiniIso"));
    }

  if( fillTauInfo_ )
    {
      if (debug__) std::cout<< " fillTauInfo_"<<std::endl;
      tauTree_                               = new hpstauInfo("",tree_, false);
      tauTree_->tauToken                     = consumes<pat::TauCollection>(iConfig.getUntrackedParameter<edm::InputTag> ("tauLabel"));
      tauTree_->theBeamSpotToken             = consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"));
    }

  if( fillPhotInfo_)
    {
      if (debug__) std::cout<< " fillPhotInfo_"<<std::endl;
      photonTree_                                 = new photonTree("", tree_);
      photonTree_->photonToken                    = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag> ("photonLabel"));
      /*
      photonTree_->phoLooseIdMapToken             = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("phoLooseIdMap"));
      photonTree_->phoMediumIdMapToken            = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("phoMediumIdMap"));
      photonTree_->phoTightIdMapToken             = consumes<edm::ValueMap<bool>>(iConfig.getParameter<edm::InputTag>("phoTightIdMap"));
      photonTree_->phoMVAValuesMapToken           = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoMVAValuesMapToken"));
      photonTree_->phoChargedIsolationToken       = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoChargedIsolationToken"));
      photonTree_->phoNeutralHadronIsolationToken = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolationToken"));
      photonTree_->phoPhotonIsolationToken        = consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("phoPhotonIsolationToken"));
      */
    }

  if( fillJetInfo_ )
    {
      if (debug__) std::cout<< "fillJetInfo_ "<<std::endl;

      std::string desc             = "THIN";
      THINjetTree_                 = new jetTree(desc,tree_,iConfig);
      THINjetTree_->jetToken       = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJets",desc.data())));
      THINjetTree_->vertexToken    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      THINjetTree_->rhoForJetToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
    }

  if( fillFATJetInfo_ )
    {
      if (debug__) std::cout<< " fillFATJetInfo_"<<std::endl;
      std::string desc            = "FAT";
      FATjetTree_                 = new jetTree(desc,tree_,iConfig);
      FATjetTree_->jetToken       = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJets",desc.data())));
      FATjetTree_->vertexToken    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      FATjetTree_->rhoForJetToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
      FATjetTree_->prunedMToken   = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJetsForPrunedMass",desc.data())));
    }

  if( fillAK4PuppiJetInfo_)
    {
      if (debug__) std::cout<< " fillAK4PuppiJetInfo_"<<std::endl;
      std::string desc            = "AK4Puppi";
      AK4PuppijetTree_                 = new jetTree(desc,tree_,iConfig);
      AK4PuppijetTree_->jetToken       = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJets",desc.data())));
      AK4PuppijetTree_->vertexToken    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      AK4PuppijetTree_->rhoForJetToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
    }

  if( fillAK8PuppiJetInfo_)
    {
      if (debug__) std::cout<< " fillAK8PuppiJetInfo_"<<std::endl;
      std::string desc            = "AK8Puppi";
      AK8PuppijetTree_                 = new jetTree(desc,tree_,iConfig);
      AK8PuppijetTree_->jetToken       = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJets",desc.data())));
      AK8PuppijetTree_->vertexToken    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      AK8PuppijetTree_->rhoForJetToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
    }

  if( fillCA15PuppiJetInfo_)
    {
      if (debug__) std::cout<< " fillCA15PuppiJetInfo_"<<std::endl;
      std::string desc            = "CA15Puppi";
      CA15PuppijetTree_                 = new jetTree(desc,tree_,iConfig);
      CA15PuppijetTree_->jetToken       = consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>(Form("%sJets",desc.data())));
      CA15PuppijetTree_->vertexToken    = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("pvSrc"));
      CA15PuppijetTree_->rhoForJetToken = consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"));
    }

    tree_->Branch("isSignal", &runOnSignal_, "runOnSignal_/O");
    tree_->Branch("mass_A", &mass_A, "mass_A/I");
    tree_->Branch("mass_a", &mass_a, "mass_a/I");
    tree_->Branch("mass_chi", &mass_chi, "mass_chi/I");
}


TreeMaker::~TreeMaker()
{
  //delete tree_;
}

void
TreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  using namespace edm;

  if( fillPUweightInfo_ ) puweight_      ->Fill(iEvent);
  if( fillEventInfo_ )    eventInfo_     ->Fill(iEvent);
  if( fillMetInfo_ )      patMetTree_    ->Fill(iEvent);
  if( fillTrigInfo_ )     patHltTree_    ->Fill(iEvent);
  if( fillFilterInfo_ )   patFilterTree_ ->Fill(iEvent);

  if( fillGenInfo_ )      genInfoTree_   ->Fill(iEvent);

  if( fillElecInfo_ )     patElecTree_   ->Fill(iEvent);
  if( fillPhotInfo_ )     photonTree_    ->Fill(iEvent);
  if( fillMuonInfo_ )     patMuTree_     ->Fill(iEvent);
  if( fillTauInfo_ )      tauTree_       ->Fill(iEvent, iSetup);



  if( fillFATJetInfo_ )   FATjetTree_    ->Fill(iEvent, iSetup);
  if( fillJetInfo_ )      THINjetTree_   ->Fill(iEvent, iSetup);
  if( fillAK4PuppiJetInfo_ ) AK4PuppijetTree_->Fill(iEvent, iSetup);
  if( fillAK8PuppiJetInfo_ ) AK8PuppijetTree_->Fill(iEvent, iSetup);
  if( fillCA15PuppiJetInfo_ ) CA15PuppijetTree_->Fill(iEvent, iSetup);
  tree_->Branch("isSignal", &runOnSignal_, "runOnSignal_/O");
  tree_->Branch("mass_A", &mass_A, "mass_A/I");
  tree_->Branch("mass_a", &mass_a, "mass_a/I");
  tree_->Fill();
}


void
TreeMaker::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup) {

}

void
TreeMaker::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup){

  if( fillGenInfo_ )
    genInfoTree_   ->GetRunInfo(iRun);

}

void
TreeMaker::beginJob(){
}



void
TreeMaker::endJob() {

}

void TreeMaker::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const& iSetup)
{
  edm::Handle<GenLumiInfoHeader> gen_header;
  iLumi.getByToken(genLumiHeaderToken_, gen_header);
  runOnSignal_ = runOnrp2HDM_ || runOnrpZpB_;
  if (runOnSignal_)
  {
    scanId_ = gen_header->configDescription();
    std::cout << " scanId_ = " << scanId_ << std::endl;
    //  splitting the scanId_ string to get ma_ and mA_
    mp_list.clear();
    chi_mass.clear();
    Zp_mass.clear();
    const char delim = '_';
    tokenize(scanId_, delim, mp_list);
    if (runOnrpZpB_)
    {
      //std::cout << "mp_list  "<<mp_list[2] << mp_list[3]<<std::endl;
      const char delim2 = 'p';
      tokenize(mp_list[2], delim2, Zp_mass); //UNcomment this line for ZpBaryonic
      const char delim3 = 'i';
      tokenize(mp_list[3], delim3, chi_mass);
      mass_A = std::stoi(Zp_mass[1]);
      mass_chi = std::stoi(chi_mass[1]);
      //std::cout << "mass_A"  << mass_A << "mass_chi" << mass_chi << std::endl;
    };
    if (runOnrp2HDM_)
    {
      mass_A = std::stoi(mp_list[1]);
      mass_a = std::stoi(mp_list[3]);
    };
  }
  else
  {
    mass_A = 0;
    mass_a = 0;
    mass_chi = 0;
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TreeMaker);
