#include "ExoPieElement/TreeMaker/interface/EventCounter.h"


EventCounter::EventCounter (const edm::ParameterSet& pset){

  edm::Service<TFileService> fs;
  totalEvents_ = fs->make<TH1F>("totalEvents","totalEvents",1,0., 1.);
}

EventCounter::~EventCounter(){}


bool EventCounter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup){
  totalEvents_->Fill(0.5);
  return true;
}

DEFINE_FWK_MODULE(EventCounter);
