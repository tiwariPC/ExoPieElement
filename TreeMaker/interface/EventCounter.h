#ifndef EventCounter_h
#define EventCounter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1F.h"




class EventCounter : public edm::EDFilter{
  
 public:
  
  //! ctor
  EventCounter (const edm::ParameterSet& pset);
  
  //! dtor 
  ~EventCounter();
  
  bool filter(edm::Event& iEvent, const edm::EventSetup& iSetup);
  
  
 private:
  
  TH1F* totalEvents_;
  
};

#endif
