// Updated By : Raman Khurana 
// Dated      : Mon May 25 15:40:47 CDT 2015
// Added possible triggers for DM analysis, Jets and MET
#include "ExoPieElement/TreeMaker/interface/patFilters.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h" 
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"

patFilters::patFilters(std::string name,TTree* tree,const edm::ParameterSet& iConfig):
  baseTree(name,tree),
  runOn2016_(iConfig.getParameter<bool>("runOn2016")),
  nfilters_(0)
{
  SetBranches();
}

void
patFilters::Fill(const edm::Event& iEvent)
{
  Clear();
  using namespace edm;

  bool    _passecalBadCalibFilterUpdate = true;//this is for 2017,18
  if (!runOn2016_){
  edm::Handle< bool > passecalBadCalibFilterUpdate ;
  iEvent.getByToken(ecalBadCalibFilterUpdate_token,passecalBadCalibFilterUpdate);
  bool    _passecalBadCalibFilterUpdate =  (*passecalBadCalibFilterUpdate );  
  }

  edm::Handle<edm::TriggerResults> trigResults;
  if (not iEvent.getByToken(filterTrigResultsToken, trigResults)) {
    std::cout << ">>> TRIGGER collection for filters does not exist !!!\n";
    return;
  }
  
  
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);
  //const std::vector<std::string> & triggerNames_ = trigNames.triggerNames();
  //std::cout <<" filter size = "<<trigResults->size()<<std::endl;
  for (unsigned int i=0; i<trigResults->size(); i++)
    {
      std::string trigName = trigNames.triggerName(i);
      // lepton triggers
      size_t foundallFlag=trigName.find("Flag_");
            

      if(false) std::cout<<" trigName = "<<trigName
			<<" : "<<trigResults->accept(i)
			<<std::endl;
      
      if ( foundallFlag==std::string::npos  )     	continue;
      
            
      filterName_.push_back(trigName);
      bool trigResult = trigResults->accept(i); //bool not to use
      filterResult_.push_back(trigResult);
      nfilters_++;
    }
  if (!runOn2016_){
  filterName_.push_back("ecalBadCalibReducedMINIAODFilter");
  filterResult_.push_back(_passecalBadCalibFilterUpdate);
  }
}

void patFilters::SetBranches(){
  
  AddBranch(&nfilters_,"nfilters");
  AddBranch(&filterResult_,"filterResult");
  AddBranch(&filterName_,"filterName");
  
}

void
patFilters::Clear(){
  nfilters_ = 0;
  filterResult_.clear();
  filterName_.clear();
}


