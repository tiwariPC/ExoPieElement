#include "ExoPieElement/TreeMaker/interface/patMetTree.h"
patMetTree::patMetTree(std::string name, TTree* tree):
  baseTree(name,tree)
{
  SetBranches();
}


patMetTree::~patMetTree(){
}


void
patMetTree::Fill(const edm::Event& iEvent){
  Clear();

  // adding Raw PF MET to the tree
  edm::Handle<reco::PFMETCollection> patMetRawHandle;
  if(not iEvent.getByToken(pfMETRawToken,patMetRawHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"pfMetRaw" <<std::endl; exit(0);}

  // adding Type-1 MET to the tree
  edm::Handle<pat::METCollection> patMetHandle;
  if(not iEvent.getByToken(pfMETToken,patMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"pfMet"<<std::endl; exit(0);}

  // adding PUPPI MET
  edm::Handle<pat::METCollection> puppiMetHandle;
  if(not iEvent.getByToken(puppimetToken, puppiMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
             <<"slimmedMETsPuppi"<<std::endl; exit(0);}
  //slimmedMETsPuppi

  auto metraw=patMetRawHandle.product()->begin();
  patMetRawPt_ = metraw->et();
  patMetRawPhi_ = metraw->phi();
  patMetRawSumEt_ = metraw->sumEt();
  patMetRawCov00_ = metraw->getSignificanceMatrix()(0,0);
  patMetRawCov01_ = metraw->getSignificanceMatrix()(0,1);
  patMetRawCov10_ = metraw->getSignificanceMatrix()(1,0);
  patMetRawCov11_ = metraw->getSignificanceMatrix()(1,1);

  pat::METCollection::const_iterator met=patMetHandle.product()->begin();
  patMetCorrPt_    = met->et();
  patMetCorrPhi_   = met->phi();
  patMetCorrSumEt_ = met->sumEt();
  patMetCorrSig_   = met->significance() < 1.e10 ? met->significance() : 0;


  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::MuonEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::MuonEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::ElectronEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::ElectronEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::TauEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::TauEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::NoShift));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::METUncertaintySize));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResUpSmear));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResDownSmear));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::METFullUncertaintySize));


  auto metpuppi = puppiMetHandle.product()->begin();
  puppiMETPt_         = metpuppi->et();
  puppiMETPhi_        = metpuppi->phi();
  puppiMETSumEt_      = metpuppi->sumEt();
  puppiMETSig_        = metpuppi->significance() < 1.e10 ? met->significance() : 0 ;

  /*  JetResUp=0, JetResDown=1, JetEnUp=2, JetEnDown=3,
    MuonEnUp=4, MuonEnDown=5, ElectronEnUp=6, ElectronEnDown=7,
    TauEnUp=8, TauEnDown=9, UnclusteredEnUp=10, UnclusteredEnDown=11,
    PhotonEnUp=12, PhotonEnDown=13, NoShift=14, METUncertaintySize=15,
    JetResUpSmear=16, JetResDownSmear=17, METFullUncertaintySize=18
  */
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetEnDown));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::MuonEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::MuonEnDown));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::ElectronEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::ElectronEnDown));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::TauEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::TauEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::NoShift));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::METUncertaintySize));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResUpSmear));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResDownSmear));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::METFullUncertaintySize));

}
bool met_extra = false;
void
patMetTree::SetBranches(){

  AddBranch(&patMetCorrPt_, "MetCorrPt");
  AddBranch(&patMetCorrPhi_, "MetCorrPhi");
  AddBranch(&patMetCorrUnc_, "MetCorrUnc");
  if (met_extra){
    AddBranch(&patMetCorrSumEt_, "MetCorrSumEt");
    AddBranch(&patMetCorrSig_, "MetCorrSig");

    AddBranch(&patMetRawPt_, "MetRawPt");
    AddBranch(&patMetRawPhi_, "MetRawPhi");
    AddBranch(&patMetRawSumEt_, "MetRawSumEt");
    AddBranch(&patMetRawCov00_, "MetRawCov00");
    AddBranch(&patMetRawCov01_, "MetRawCov01");
    AddBranch(&patMetRawCov10_, "MetRawCov10");
    AddBranch(&patMetRawCov11_, "MetRawCov11");

    AddBranch(&mvaMetPt_,     "mvaMetPt");
    AddBranch(&mvaMetPhi_,    "mvaMetPhi");
    AddBranch(&mvaMetSumEt_,  "mvaMetSumEt");
    AddBranch(&mvaMetSig_,    "mvaMetSig");

    AddBranch(&puppiMETPt_,     "puppiMETPt");
    AddBranch(&puppiMETPhi_,    "puppiMETPhi");
    AddBranch(&puppiMETSumEt_,  "puppiMETSumEt");
    AddBranch(&puppiMETSig_,    "puppiMETSig");
    AddBranch(&puppiMETUnc_,    "puppiMETUnc");
  }
}


void
patMetTree::Clear(){

  float dummy = -99999;
  patMetCorrPt_= dummy;
  patMetCorrPhi_= dummy;
  patMetCorrSumEt_= dummy;
  patMetCorrSig_= dummy;
  patMetCorrUnc_.clear();

  patMetRawPt_= dummy;
  patMetRawPhi_= dummy;
  patMetRawSumEt_= dummy;
  patMetRawCov00_= dummy;
  patMetRawCov01_= dummy;
  patMetRawCov10_= dummy;
  patMetRawCov11_= dummy;

  mvaMetPt_     = dummy ;
  mvaMetPhi_    = dummy ;
  mvaMetSumEt_  = dummy ;
  mvaMetSig_    = dummy ;


  puppiMETPt_     = dummy ;
  puppiMETPhi_    = dummy ;
  puppiMETSumEt_  = dummy ;
  puppiMETSig_    = dummy ;
  puppiMETUnc_.clear();



}
