#ifndef __ELEC_TREE_H_
#define __ELEC_TREE_H_
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TClonesArray.h"
#include <bitset>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"
#include "DataFormats/Common/interface/RefHolder.h"
#include "DataFormats/Common/interface/RefVectorHolder.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"


#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"


#include "ExoPieElement/TreeMaker/interface/baseTree.h"
#include "ExoPieElement/TreeMaker/interface/utils.h"




#include "TTree.h"
#include "Math/VectorUtil.h"

using namespace std;
using namespace edm;
class patElecTree : public baseTree {
 public:
  patElecTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patElecTree();
  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<reco::VertexCollection>          vertexToken;
  edm::EDGetTokenT<double>                          rhoForLepToken;
  edm::EDGetTokenT<edm::View<pat::Electron>>        eleToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection>  pfCandToken;

  /*
  // cut-based
  edm::EDGetTokenT<edm::ValueMap<bool>> eleVetoIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleLooseIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMediumIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleHEEPIdMapToken;

  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> eleVetoIdCFToken;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> eleLooseIdCFToken;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> eleMediumIdCFToken;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> eleTightIdCFToken;
  edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>> eleHEEPIdCFToken;
  

  // MVA based
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMVAMediumIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMVATightIdMapToken;
  
  // MVA values and categories (optional)
  edm::EDGetTokenT<edm::ValueMap<float>> mvaValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<int>> mvaCategoriesMapToken;
  */

 private:

  //TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patElecTree();
  void SetBranches();

  double r_iso_min_;
  double r_iso_max_;
  double kt_scale_;
  bool charged_only_;
  EffectiveAreas eAreasElectrons;

  // ntuple variables 
  float patElecRho_;
  int nEle_;

  TClonesArray *patElecP4_;

  std::vector<float> patElecPx_;
  std::vector<float> patElecPy_;
  std::vector<float> patElecPz_;
  std::vector<float> patElecE_;


  std::vector<bool> patElecInBarrel_;
  std::vector<bool> patElecInEndcap_;

  std::vector<int> patElecCharge_;
  std::vector<int> patElecChargeConsistent_;

  std::vector<float> patElecaloEnergy_;

  std::vector<float> patElecScEt_;
  std::vector<float> patElecScEn_;
  std::vector<float> patElecScPreEn_;
  std::vector<float> patElecScEta_;
  std::vector<float> patElecScPhi_;
  std::vector<float> patElecScRawEn_;
  std::vector<float> patElecScEtaWidth_;
  std::vector<float> patElecScPhiWidth_;

  std::vector<float> patElecR9_;
  std::vector<float> patElecHoverE_;

  std::vector<float> patElecD0_;
  std::vector<float> patElecDz_;

  std::vector<float> patElecEoverP_;
  std::vector<float> patElecBrem_;
  std::vector<float> patElecdEtaAtVtx_;
  std::vector<float> patElecdPhiAtVtx_;
  std::vector<float> patElecSigmaIEtaIEta_;
  std::vector<float> patElecSigmaIEtaIPhi_;
  std::vector<float> patElecSigmaIPhiIPhi_;

  std::vector<bool> patElecConvVeto_;
  std::vector<int> patElecMissHits_;
  std::vector<float> patElecEoverPInv_;

  std::vector<float> patElecdEtaseedAtVtx_;
  std::vector<float> patElecE1x5_;
  std::vector<float> patElecE2x5_;
  std::vector<float> patElecE5x5_;

  std::vector<float> patElecSigmaIEtaIEtaFull5x5_;
  std::vector<float> patElecE1x5Full5x5_;
  std::vector<float> patElecE2x5Full5x5_;
  std::vector<float> patElecE5x5Full5x5_;
  std::vector<float> patElecR9Full5x5_;

  std::vector<float> patElecChHadIso_;
  std::vector<float> patElecNeHadIso_;
  std::vector<float> patElecGamIso_;
  std::vector<float> patElecPUPt_;

  // for MVA preselection
  std::vector<float> patElecEcalPFClusterIso_;
  std::vector<float> patElecHcalPFClusterIso_;  

  // miniIso input
  std::vector<float> patElecMiniIso_ch_;
  std::vector<float> patElecMiniIso_nh_;
  std::vector<float> patElecMiniIso_ph_;
  std::vector<float> patElecMiniIso_pu_;
  std::vector<float> patElecMiniIso_r_;
  std::vector<float> patElecMiniIsoBeta_; // subtracting pu via beta
  std::vector<float> patElecMiniIsoEA_; // subtracting pu via rho/EA


  std::vector<bool> patElecEcalDrivenSeed_;
  std::vector<bool> patElecEcalDriven_;
  std::vector<float> patElecDr03EcalRecHitSumEt_;
  std::vector<float> patElecDr03HcalDepth1TowerSumEt_;
  std::vector<float> patElecDr03HcalDepth2TowerSumEt_;
  std::vector<float> patElecDr03HcalTowerSumEt_;
  std::vector<float> patElecDr03TkSumPt_;



  std::vector<bool> isPassVeto_;
  std::vector<bool> isPassLoose_;
  std::vector<bool> isPassMedium_;
  std::vector<bool> isPassTight_;
  std::vector<bool> isPassHEEP_;
  std::vector<bool> isPassVetoNoIso_;
  std::vector<bool> isPassLooseNoIso_;
  std::vector<bool> isPassMediumNoIso_;
  std::vector<bool> isPassTightNoIso_;
  std::vector<bool> isPassHEEPNoIso_;
  std::vector<bool> isPassMVAMedium_;
  std::vector<bool> isPassMVATight_;
  

  std::vector<float> mvaValue_;
  std::vector<int>   mvaCategory_;
  
  


};
#endif
