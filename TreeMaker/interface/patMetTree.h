#ifndef __MET_TREE_H_
#define __MET_TREE_H_
/*
  Updated By: Raman Khurana 
  Date      : 24 June 2015.

 -- it can save three types of MET now : 
   = PFMET uncorrected
   = PFMET corrected
   = MVA MET

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
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

#include "ExoPieElement/TreeMaker/interface/utils.h"
#include "ExoPieElement/TreeMaker/interface/baseTree.h"
using namespace std;
using namespace edm;



class patMetTree : public baseTree{

 public:
  patMetTree(std::string name, TTree* tree);
  ~patMetTree();
  void Fill(const edm::Event& iEvent);
  void Clear();

  edm::EDGetTokenT<reco::PFMETCollection>           pfMETRawToken;
  edm::EDGetTokenT<pat::METCollection>              pfMETToken;
  edm::EDGetTokenT<reco::PFMETCollection>           pfMVAMETToken;
  edm::EDGetTokenT<pat::METCollection>             puppimetToken;
  
  
 private:

  patMetTree(){};
  void SetBranches();
    
  float patMetCorrPt_;  
  float patMetCorrPhi_; 
  float patMetCorrSumEt_;
  float patMetCorrSig_;
  std::vector<float> patMetCorrUnc_;

  float patMetRawPt_;
  float patMetRawPhi_;
  float patMetRawSumEt_;
  float patMetRawCov00_;
  float patMetRawCov01_;
  float patMetRawCov10_;
  float patMetRawCov11_;


  float mvaMetPt_;
  float mvaMetPhi_;
  float mvaMetSumEt_;
  float mvaMetSig_;

  float puppiMETPt_; 
  float puppiMETPhi_; 
  float puppiMETSumEt_; 
  float puppiMETSig_; 
  std::vector<float>  puppiMETUnc_; 
  
};

#endif

