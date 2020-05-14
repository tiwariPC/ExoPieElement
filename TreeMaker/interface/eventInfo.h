#ifndef __eventInfo__
#define __eventInfo__

/*
  Updated by: Shin-Shan Yu
  Date      : 20 March 2016
  Replace getByLabel with getByToken
*/

#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "ExoPieElement/TreeMaker/interface/baseTree.h"

class eventInfo : public baseTree {

 public:
  eventInfo(std::string name, TTree* tree);
  ~eventInfo();

  void Fill(const edm::Event& iEvent);
  void Clear();
  edm::EDGetTokenT<reco::VertexCollection>          vertexToken;

 private:

  eventInfo(){};
  void SetBranches();

  bool isData_;

  unsigned long nEvt_;
  unsigned long nRun_;
  unsigned long nLumiS_;
  unsigned long bunchX_;

  int nVtx_;

  TClonesArray *vertexP3_;
  float _prefiringweight;
  float _prefiringweightup;
  float _prefiringweightdown;
};

#endif
