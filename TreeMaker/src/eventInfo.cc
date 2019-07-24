#include <iostream>
#include "ExoPieElement/TreeMaker/interface/eventInfo.h"
#include "TVector3.h"

eventInfo::eventInfo(std::string name, TTree* tree):
  baseTree(name, tree)
{
  vertexP3_ =   new TClonesArray("TVector3");
  SetBranches();

}


eventInfo::~eventInfo(){
  delete vertexP3_;
}


void
eventInfo::Fill(const edm::Event& iEvent){
  Clear();

  isData_ = iEvent.isRealData();
  nEvt_   = iEvent.id().event();
  nRun_   = iEvent.id().run();
  nLumiS_ = iEvent.luminosityBlock();
  bunchX_ = iEvent.bunchCrossing();


  edm::Handle<reco::VertexCollection> recVtxs_;

  if (iEvent.getByToken(vertexToken, recVtxs_)) {
    for (size_t i=0; i<recVtxs_->size(); ++i) {
      
      //bool isFake = ((*recVtxs_)[i].chi2()==0 && (*recVtxs_)[i].ndof()==0);
      if((*recVtxs_)[i].ndof() >= 4 && fabs((*recVtxs_)[i].z()) <= 24 && fabs((*recVtxs_)[i].position().rho()) <= 2
	 //&& !(isFake)
	 )
	{
	  nVtx_++;
	  TVector3 v3((*recVtxs_)[i].x(),
		      (*recVtxs_)[i].y(),
		      (*recVtxs_)[i].z());
	  new( (*vertexP3_)[nVtx_-1]) TVector3(v3);
	} // if satifying good vertices
    }
  }
}

void
eventInfo::SetBranches(){
  AddBranch(&isData_, "isData");
  AddBranch(&nEvt_,"eventId");
  AddBranch(&nRun_,  "runId");
  AddBranch(&nLumiS_, "lumiSection");
  AddBranch(&bunchX_, "bunchXing");
  AddBranch(&nVtx_, "nVtx");
  AddBranch(&vertexP3_,"vertexP3");
}


void 
eventInfo::Clear(){

  isData_  = false;
  nEvt_   = 0;
  nRun_   = 0;
  nLumiS_ = 0;
  bunchX_ = 0;
  nVtx_ = 0;
  vertexP3_->Clear();
}

