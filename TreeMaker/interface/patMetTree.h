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

  edm::EDGetTokenT<pat::METCollection>              pfMETToken;
  edm::EDGetTokenT<pat::METCollection>              pfMETModifiedToken;
  edm::EDGetTokenT<pat::METCollection>             puppimetToken;


 private:

  patMetTree(){};
  void SetBranches();

  float patMetCorrPt_;
  float patMetCorrPhi_;
  float patMetCorrSumEt_;
  float patMetCorrSig_;
  std::vector<float> patMetCorrUnc_;

  float patmodifiedMetCorrPt_;
  float patmodifiedMetCorrPhi_;
  float patmodifiedMetCorrSumEt_;
  float patmodifiedMetCorrSig_;
  std::vector<float> patmodifiedMetCorrUnc_;

  float patMet_smear_;
  float patmodifiedMet_smear_;

  float patCaloMETPt_;
  float patCaloMETPhi_;
  float patCaloMETSumEt_;

  float patGenMETPt_;
  float patGenMETPhi_;
  float patGenMETSumEt_;

  float patMetRawPt_;
  float patMetRawPhi_;
  float patMetRawSumEt_;


  float puppiMETPt_;
  float puppiMETPhi_;
  float puppiMETSumEt_;
  float puppiMETSig_;


  float CHSMETPt_;
  float CHSMETPhi_;
  float CHSMETSumEt_;

  float TRKMETPt_;
  float TRKMETPhi_;
  float TRKMETPSumEt_;





  std::vector<float>  puppiMETUnc_;
  bool is_Data;

};

#endif
