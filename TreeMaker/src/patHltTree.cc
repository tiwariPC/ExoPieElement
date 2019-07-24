// Updated By : Raman Khurana, Shin-Shan Eiko Yu
// Dated      : Mon May 25 15:40:47 CDT 2015
// Added possible triggers for DM analysis, Jets and MET
#include <memory>
#include <cmath>
#include <vector>
#include <string.h>
#include <TString.h>
#include "ExoPieElement/TreeMaker/interface/patHltTree.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"



patHltTree::patHltTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(name,tree),
  saveAllTrigPaths_(iConfig.getParameter<bool>("saveAllTrigPaths")),
  nTrigs_(0),
  nTrigObj_(0)
{
    SetBranches();
}

void
patHltTree::Fill(const edm::Event& iEvent)
{
  Clear();
  using namespace edm;

  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  iEvent.getByToken(triggerPrescalesToken,triggerPrescales);

  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjectsToken,triggerObjects);

  edm::Handle<edm::TriggerResults> trigResults;
  if (not iEvent.getByToken(trigResultsToken, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return;
  }
std::vector<std::string> triggerlist={"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight","HLT_PFMETNoMu140_PFMHTNoMu140_IDTight","HLT_IsoMu27","HLT_IsoTkMu27","HLT_IsoMu24","HLT_IsoTkMu24","HLT_Ele27_WPTight_Gsf","HLT_Ele35_WPTight_Gsf","HLT_Photon200", "HLT_Ele115_CaloIdVT_GsfTrkIdT_v","HLT_Ele32_WPTight_Gsf_v"};
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);

  for (unsigned int i=0; i<trigResults->size(); i++)
    {
      std::string trigName = trigNames.triggerName(i);


      bool trigResult = trigResults->accept(i); //bool not to use
      // if(!trigResult && !saveAllTrigPaths_)continue;
      int prescale = triggerPrescales->getPrescaleForIndex(i);
      if(prescale!=1 && !saveAllTrigPaths_)continue;

      if(false) std::cout<<" trigName = "<<trigName
       			<<" : " << trigResults->accept(i)
			<<" : " << triggerPrescales->getPrescaleForIndex(i)
       			<<std::endl;
      std::string trigName_123 = trigName.substr(0, trigName.find("_v"));
      if (std::find(triggerlist.begin(), triggerlist.end(), trigName_123) == triggerlist.end()) continue;
      trigName_.push_back(trigName);
      trigResult_.push_back(trigResult);
      trigPrescale_.push_back(prescale);
      nTrigs_++;
    }
  //std::cout<<"nTrigs: "<< nTrigs_<<std::endl;
  for (pat::TriggerObjectStandAlone trigObject : *triggerObjects) {
      trigObject.unpackFilterLabels(iEvent, *trigResults);
      trigObject.unpackPathNames(trigNames);
      for (unsigned h = 0; h < trigObject.filterLabels().size(); ++h) {
          //std::cout << " " << trigObject.filterLabels()[h];
          std::string fullpath1 = trigObject.filterLabels()[h];
          std::string fullpath2 = std::to_string(nTrigObj_);
          std::string fullpath = fullpath1 + "_Obj_" + fullpath2;
          trigFilterLabels_.push_back(fullpath);
      }
      std::vector<std::string> pathNamesAll = trigObject.pathNames(false);
      std::vector<std::string>  pathNamesLast = trigObject.pathNames(true);

      //std::cout << "\t   Paths (" << pathNamesAll.size()<<"/"<<pathNamesLast.size()<<"):    ";
      for (unsigned h = 0, n = pathNamesAll.size(); h < n; ++h) {
          std::string fullpathName1 = pathNamesAll[h];
          std::string fullpathName2 = std::to_string(nTrigObj_);
          std::string fullpathName = fullpathName1 + "_Obj_" + fullpathName2;
          trigPathName_.push_back(fullpathName);
      }

      trigObj_pT_.push_back(trigObject.pt());
      trigObj_eta_.push_back(trigObject.eta());
      trigObj_phi_.push_back(trigObject.phi());
      nTrigObj_++;
  }
  //std::cout<<"# trigObjects: "<<nTrigObj_<<std::endl;

}

void patHltTree::SetBranches(){

  AddBranch(&nTrigs_,"nTrigs");
  AddBranch(&nTrigObj_,"nTrigObj");
  AddBranch(&trigResult_,"trigResult");
  AddBranch(&trigName_,"trigName");
  AddBranch(&trigPrescale_,"trigPrescale");
  AddBranch(&trigObj_pT_,"trigObj_pT");
  AddBranch(&trigObj_eta_,"trigObj_eta");
  AddBranch(&trigObj_phi_,"trigObj_phi");
  AddBranch(&trigFilterLabels_,"trigFilterLabels");
  AddBranch(&trigPathName_,"trigPathName");

}

void
patHltTree::Clear(){
  nTrigs_ = 0;
  nTrigObj_ = 0;
  trigResult_.clear();
  trigName_.clear();
  trigPrescale_.clear();
  trigObj_pT_.clear();
  trigObj_eta_.clear();
  trigObj_phi_.clear();
  trigFilterLabels_.clear();
  trigPathName_.clear();
}
