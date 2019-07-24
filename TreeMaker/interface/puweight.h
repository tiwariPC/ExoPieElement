#ifndef __puweight__
#define __puweight__

/*
  Updated by: Shin-Shan Yu
  Date      : 20 March 2016
  Replace getByLabel with getByToken
*/

#include <memory>
#include <string>
#include <iostream>
#include "TTree.h" 
#include "FWCore/Framework/interface/Event.h" 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "ExoPieElement/TreeMaker/interface/baseTree.h"

class puweight : public baseTree {

 public:

  puweight(std::string name, TTree* tree);
  ~puweight();

  void Fill(const edm::Event& iEvent); 
  void Clear();

  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken;

 private:

  puweight(){};
  void SetBranches();


  float nTrueInt_;
  int   nPUVert_;

};

#endif

