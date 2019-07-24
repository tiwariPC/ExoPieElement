#ifndef __MUON_TREE_H_
#define __MUON_TREE_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>

#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"


#include "ExoPieElement/TreeMaker/interface/baseTree.h"
#include "ExoPieElement/TreeMaker/interface/utils.h"

using namespace std;
using namespace edm;

class patMuonTree : public baseTree {

 public:

  patMuonTree(std::string name, TTree* tree, const edm::ParameterSet& cfg);
  ~patMuonTree();

  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<reco::VertexCollection>          vertexToken;
  edm::EDGetTokenT<double>                          rhoForLepToken;
  edm::EDGetTokenT<pat::MuonCollection>             muToken;
  edm::EDGetTokenT<pat::PackedCandidateCollection>  pfCandToken;


 private:
  //TTree* tree_;
  //Dont Allow User to Call the Default Constructor.
  patMuonTree();
  void SetBranches();

  double r_iso_min_;
  double r_iso_max_;
  double kt_scale_;
  bool charged_only_;
  EffectiveAreas eAreasMuons;

  // ntuple variabes
  int nMu;
  TClonesArray *patMuonP4;

  std::vector<float> patMuonPx_;
  std::vector<float> patMuonPy_;
  std::vector<float> patMuonPz_;
  std::vector<float> patMuonE_;


  std::vector<int> patMuonType;
  std::vector<int> patMuonCharge;


  std::vector<bool> isGlobalMuon;
  std::vector<bool> isTrackerMuon;
  std::vector<bool> isPFMuon;
  
  std::vector<bool> isTightMuon;
  std::vector<bool> isLooseMuon;
  std::vector<bool> isMediumMuon;
  std::vector<bool> isSoftMuon;
  std::vector<bool> isHighPtMuon;
  std::vector<bool> isCustomTrackerMuon;
  
  std::vector<int>   patMuonITrkIndex;
  std::vector<int>   patMuonSegIndex;
  std::vector<int>   patMuonNSeg;
  std::vector<bool>   patMuonGood;
  std::vector<bool>   patMuonIsGood;

  std::vector<float> patMuonTrkPt;
  std::vector<float> patMuonTrkPtErr;
  std::vector<float> patMuondxy;
  std::vector<float> patMuondz;
  std::vector<float> patMuonsegmentCompatibility;
  std::vector<float> patMuonchi2LocalPosition;
  std::vector<float> patMuontrkKink;

  std::vector<float> patMuonInnerdxy;
  std::vector<float> patMuonInnerdz;
  std::vector<int>   patMuonTrkLayers;
  std::vector<int>   patMuonPixelLayers;
  std::vector<int>   patMuonPixelHits;
  std::vector<int>   patMuonHits;
  std::vector<float> patMuonTrkQuality;
  std::vector<float> patMuonChi2NDF;
  std::vector<float> patMuonInnervalidFraction;
  std::vector<int>   patMuonMatches;

  std::vector<float> patMuonTrkIso;
  std::vector<float> patMuonHcalIso;
  std::vector<float> patMuonEcalIso;
  std::vector<float> patMuonChHadIso;
  std::vector<float> patMuonNeHadIso;
  std::vector<float> patMuonGamIso;
  std::vector<float> patMuonPUPt;
  std::vector<float> patMuonInnerTrkPt;

  // miniIso input
  std::vector<float> patMuonMiniIso_ch;
  std::vector<float> patMuonMiniIso_nh;
  std::vector<float> patMuonMiniIso_ph;
  std::vector<float> patMuonMiniIso_pu;
  std::vector<float> patMuonMiniIso_r;
  std::vector<float> patMuonMiniIsoBeta; // subtracting pu via beta
  std::vector<float> patMuonMiniIsoEA; // subtracting pu via rho/EA

 
  
};
#endif
