#ifndef patHlttree_h
#define patHlttree_h

/*
  Updated by: Shin-Shan Yu
  Date      : 20 March 2016
  Replace getByLabel with getByToken
*/

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "ExoPieElement/TreeMaker/interface/baseTree.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

class patHltTree : public baseTree{

 public:
  patHltTree(std::string name,TTree* tree,const edm::ParameterSet& iConfig);
  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<edm::TriggerResults>             trigResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>     triggerPrescalesToken;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>     triggerObjectsToken;

 private:

  patHltTree(){};
  void SetBranches();
  bool runOn2017_;
  bool runOn2016_;
  bool saveAllTrigPaths_;
  int nTrigs_;
  int nTrigObj_;
  std::vector<bool> trigResult_;
  std::vector<std::string> trigName_;
  std::vector<int> trigPrescale_;
  std::vector<std::string> trigFilterLabels_;
  std::vector<std::string> trigPathName_;
  std::vector<float> trigObj_pT_;
  std::vector<float> trigObj_eta_;
  std::vector<float> trigObj_phi_;
  edm::InputTag trigTag;
};

#endif
