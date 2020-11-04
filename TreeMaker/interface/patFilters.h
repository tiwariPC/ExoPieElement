#ifndef patFilters_h
#define patFilters_h

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
#include "ExoPieElement/TreeMaker/interface/baseTree.h"

class patFilters : public baseTree{

 public:
  patFilters(std::string name,TTree* tree,const edm::ParameterSet& iConfig);
  void Fill(const edm::Event& iEvent);
  void Clear();
  edm::EDGetTokenT<bool>                            HBHETToken;
  edm::EDGetTokenT<bool>                            HBHELToken;
  edm::EDGetTokenT<bool>                            HBHEIsoToken;
  edm::EDGetTokenT<edm::TriggerResults>             filterTrigResultsToken;
  edm::EDGetTokenT< bool >ecalBadCalibFilterUpdate_token ;
  edm::EDGetTokenT<bool> BadChCandFilterToken_;
  edm::EDGetTokenT<bool> BadPFMuonFilterToken_;
  edm::EDGetTokenT<bool> BadGlobalMuonFilterToken_;
  edm::EDGetTokenT<bool> CloneGlobalMuonFilterToken_;

 private:
  
  patFilters(){};
  void SetBranches();

  bool runOn2016_;
  int nfilters_;
  std::vector<bool> filterResult_;
  bool hbhet_;
  bool hbhel_;
  bool hbheIso_;
  bool  filterbadChCandidate;
  bool filterbadPFMuon;
  bool filterbadGlobalMuon;
  bool filtercloneGlobalMuon;
  std::vector<std::string> filterName_;
  edm::InputTag filterTag;
};

#endif
