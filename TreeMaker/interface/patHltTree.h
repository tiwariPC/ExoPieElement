#ifndef patHlttree_h
#define patHlttree_h

/*
  Updated by: Shin-Shan Yu
  Date      : 20 March 2016
  Replace getByLabel with getByToken
*/

#include<iostream>
#include<string>
#include<vector>

#include "TTree.h"
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "ExoPieElement/TreeMaker/interface/baseTree.h"

class patHltTree : public baseTree{

 public:
  patHltTree(std::string name,TTree* tree,const edm::ParameterSet& iConfig);
  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<edm::TriggerResults>             trigResultsToken;
  edm::EDGetTokenT<pat::PackedTriggerPrescales>     triggerPrescalesToken;
  
 private:

  patHltTree(){};
  void SetBranches();
  bool saveAllTrigPaths_;
  
  int nTrigs_;
  std::vector<bool> trigResult_;
  std::vector<std::string> trigName_;
  std::vector<int> trigPrescale_;
  edm::InputTag trigTag;
};

#endif
