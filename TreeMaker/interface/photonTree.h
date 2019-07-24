#ifndef __PHOTON_TREE_H_
#define __PHOTON_TREE_H_

/*
Log:
Sep 10, 2011
Anil Singh: Empty template created. 
-----
30 December 2015
Raman Khurana: Added Photon ID Variables
-----
20 March 2016
Shin-Shan Yu: replace getByLabel with getByToken
*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "ExoPieElement/TreeMaker/interface/utils.h"
#include "ExoPieElement/TreeMaker/interface/baseTree.h"

using namespace std;
using namespace edm;

class photonTree : public baseTree{

 public:
  photonTree(std::string name, TTree* tree);
  ~photonTree();
  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<edm::View<pat::Photon>> photonToken;
  edm::EDGetTokenT<edm::ValueMap<bool>>    phoLooseIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>>    phoMediumIdMapToken;
  edm::EDGetTokenT<edm::ValueMap<bool>>    phoTightIdMapToken;

  edm::EDGetTokenT<edm::ValueMap<float>>   phoMVAValuesMapToken;
  edm::EDGetTokenT<edm::ValueMap<float>>   phoChargedIsolationToken; 
  edm::EDGetTokenT<edm::ValueMap<float>>   phoNeutralHadronIsolationToken; 
  edm::EDGetTokenT<edm::ValueMap<float>>   phoPhotonIsolationToken; 

 private:

  photonTree(){};
  void SetBranches();
  bool usePFObjects_;


  
  //variables which would become branches
  int nPho_;
  TClonesArray *photonP4_;
  std::vector<float> phoPx_;
  std::vector<float> phoPy_;
  std::vector<float> phoPz_;
  std::vector<float> phoE_;

  std::vector<bool> isPassLoose;
  std::vector<bool> isPassMedium;
  std::vector<bool> isPassTight;
  std::vector<float> phoIDMVA_;
  
  vector<float>  phoSCE_;
  vector<float>  phoSCRawE_;
  vector<float>  phoSCEta_;
  vector<float>  phoSCPhi_;
  vector<float>  phoSCEtaWidth_;
  vector<float>  phoSCPhiWidth_;
  vector<float>  phoSCBrem_;
  vector<int>    phohasPixelSeed_;
  vector<int>    phoEleVeto_;
  vector<float>  phoR9_;
  vector<float>  phoHoverE_;
  vector<float>  phoSigmaIEtaIEta_;
  vector<float>  phoSigmaIEtaIPhi_;
  vector<float>  phoSigmaIPhiIPhi_;
  vector<float>  phoSigmaIEtaIEtaFull5x5_;
  vector<float>  phoR9Full5x5_;

  vector<float>  phoPFChIso_;
  vector<float>  phoPFPhoIso_;
  vector<float>  phoPFNeuIso_;

};

#endif

