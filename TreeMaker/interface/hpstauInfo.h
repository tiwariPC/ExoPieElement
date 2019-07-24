#ifndef __HPSTAU_INFO_H_
#define __HPSTAU_INFO_H_
//#include "MVAElectronID.h"


/*
  Updated by: Shin-Shan Yu
  Date      : 20 March 2016
  Replace getByLabel with getByToken
*/

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "utils.h"
#include "baseTree.h"
#include "DataFormats/PatCandidates/interface/Tau.h"


class hpstauInfo : public baseTree{

 public:
  hpstauInfo(std::string name, TTree* tree, bool debug);
  ~hpstauInfo();
  void Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  void Clear();
  bool debug_;

  //tau tokens
  edm::EDGetTokenT<pat::TauCollection>              tauToken;
  edm::EDGetTokenT<reco::BeamSpot>                  theBeamSpotToken;


 private:

  hpstauInfo(){};
  void SetBranches();
  //variables which would become branches

  int HPSTau_n;
  TClonesArray *HPSTau_4Momentum, *HPSTau_Vposition ;
  std::vector<float> taupt;

  std::vector<float> TauPx_;
  std::vector<float> TauPy_;
  std::vector<float> TauPz_;
  std::vector<float> TauE_;


  std::vector<bool>  HPSTau_leadPFChargedHadrCand;
  std::vector<bool> HPSTau_leadPFChargedHadrCand_trackRef;

  std::vector<bool> disc_againstElectronLoose;
  std::vector<bool> disc_againstElectronMedium;
  std::vector<bool> disc_againstElectronTight;
  std::vector<bool> disc_againstElectronLooseMVA6;
  std::vector<bool> disc_againstElectronMediumMVA6;
  std::vector<bool> disc_againstElectronTightMVA6;
  std::vector<bool> disc_againstElectronVLooseMVA6;
  std::vector<bool> disc_againstElectronVTightMVA6;
  std::vector<bool> disc_againstMuonLoose;
  std::vector<bool> disc_againstMuonMedium;
  std::vector<bool> disc_againstMuonTight;
  std::vector<bool> disc_againstMuonLoose2;
  std::vector<bool> disc_againstMuonMedium2;
  std::vector<bool> disc_againstMuonTight2;
  std::vector<bool> disc_againstMuonLooseMVA;
  std::vector<bool> disc_againstMuonMediumMVA;
  std::vector<bool> disc_againstMuonTightMVA;
  std::vector<bool> disc_againstMuonLoose3;
  std::vector<bool> disc_againstMuonTight3;

  std::vector<bool> disc_byVLooseCombinedIsolationDeltaBetaCorr;
  std::vector<bool> disc_byLooseCombinedIsolationDeltaBetaCorr;
  std::vector<bool> disc_byMediumCombinedIsolationDeltaBetaCorr;
  std::vector<bool> disc_byTightCombinedIsolationDeltaBetaCorr;


  std::vector<bool> disc_byLooseIsolation;

  std::vector<bool> disc_byVLooseIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool> disc_byLooseIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool> disc_byMediumIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool> disc_byTightIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool> disc_byVTightIsolationMVArun2v1DBnewDMwLT;
  std::vector<bool> disc_byVVTightIsolationMVArun2v1DBnewDMwLT;

  std::vector<bool> disc_byIsolationMVArun2017v2DBoldDMwLTraw2017;
  std::vector<bool> disc_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byVLooseIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byTightIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byVTightIsolationMVArun2017v2DBoldDMwLT2017;
  std::vector<bool> disc_byVVTightIsolationMVArun2017v2DBoldDMwLT2017;
  // 2017 v2 new DM
  std::vector<bool> disc_byIsolationMVArun2017v2DBnewDMwLTraw2017;
  std::vector<bool> disc_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byVLooseIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byLooseIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byMediumIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byTightIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byVTightIsolationMVArun2017v2DBnewDMwLT2017;
  std::vector<bool> disc_byVVTightIsolationMVArun2017v2DBnewDMwLT2017;


  std::vector<bool> disc_byLooseCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<bool> disc_byMediumCombinedIsolationDeltaBetaCorr3Hits;
  std::vector<bool> disc_byTightCombinedIsolationDeltaBetaCorr3Hits;

  std::vector<bool> disc_decayModeFinding;
  std::vector<bool> disc_decayModeFindingNewDMs;

  std::vector<float> disc_chargedIsoPtSum;
  std::vector<float> disc_neutralIsoPtSum;
  std::vector<float> disc_puCorrPtSum;


  std::vector<float> HPSTau_NewVz;
  std::vector<int> HPSTau_charge;






};

#endif
