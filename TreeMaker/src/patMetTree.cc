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

  // adding Type-1 MET to the tree
  edm::Handle<pat::METCollection> patMetHandle;
  if(not iEvent.getByToken(pfMETToken,patMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"pfMet"<<std::endl; exit(0);}


  // adding modified Type-1 MET to the tree: EE Fixed
  edm::Handle<pat::METCollection> patMetModifiedHandle;
  if(not iEvent.getByToken(pfMETModifiedToken,patMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"modified pfMet"<<std::endl; exit(0);}

  // adding PUPPI MET
  edm::Handle<pat::METCollection> puppiMetHandle;
  if(not iEvent.getByToken(puppimetToken, puppiMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
             <<"slimmedMETsPuppi"<<std::endl; exit(0);}
  
  
  
  std::cout<<" Type 1 corrected MET default in miniaod"<<std::endl;
  
  pat::METCollection::const_iterator met=patMetHandle.product()->begin();
  // Type 1 corrected MET default in miniaod 
  patMetCorrPt_    = met->et();
  patMetCorrPhi_   = met->phi();
  patMetCorrSumEt_ = met->sumEt();
  patMetCorrSig_   = met->significance() < 1.e10 ? met->significance() : 0;

  


  // uncorrected raw met 
  patMetRawPt_ = met->uncorPt();
  patMetRawPhi_ = met->uncorPhi();
  patMetRawSumEt_ = met->uncorSumEt();
 
  // gen met  :: set is using mc flag :: in testing phase
  if (false){
  patGenMETPt_ = met->genMET()->et();
  patGenMETPhi_ = met->genMET()->phi();
  patGenMETSumEt_ = met->genMET()->sumEt();
  }
   
  
  // calo met 
  patCaloMETPt_ = met->caloMETPt();
  patCaloMETPhi_ = met->caloMETPhi();
  patCaloMETSumEt_ = met->caloMETSumEt();


  // met uncertainties, need to be changed when time comes, right now don't have much time to edit this part
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::NoShift));

  std::cout<<" after t1 met unc "<<std::endl;
  
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::MuonEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::MuonEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::ElectronEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::ElectronEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::TauEnUp));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::TauEnDown));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::METUncertaintySize));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResUpSmear));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResDownSmear));
  //patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::METFullUncertaintySize));
  
  


  // Modified Type 1 corrected MET default in miniaod :: Needed only for 2017 data mc. 
  /*
  pat::METCollection::const_iterator metmodified=patMetModifiedHandle.product()->begin();
  patMetCorrPt_    = metmodified->et();
  patMetCorrPhi_   = metmodified->phi();
  patMetCorrSumEt_ = metmodified->sumEt();
  */
  
  std::cout<< "Modified Type 1 corrected MET default in miniaod :: Needed only for 2017 data mc"<<std::endl;
// puppi met,  present in miniaod 

  auto metpuppi = puppiMetHandle.product()->begin();
  puppiMETPt_         = metpuppi->et();
  puppiMETPhi_        = metpuppi->phi();
  puppiMETSumEt_      = metpuppi->sumEt();
  puppiMETSig_        = metpuppi->significance() < 1.e10 ? met->significance() : 0 ;
  
  
  std::cout<<" puppi met,  present in miniaod "<<std::endl;
  
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
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::NoShift));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::MuonEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::MuonEnDown));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::ElectronEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::ElectronEnDown));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::TauEnUp));
  //puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::TauEnDown));
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
  AddBranch(&patMetCorrSumEt_, "MetCorrSumEt");
  
  AddBranch(&patMetCorrUnc_, "MetCorrUnc");
  AddBranch(&patMetCorrSig_, "MetCorrSig");
  
  AddBranch(&patMetRawPt_, "MetRawPt");
  AddBranch(&patMetRawPhi_, "MetRawPhi");
  AddBranch(&patMetRawSumEt_, "MetRawSumEt");

  
  AddBranch(&patGenMETPt_,     "patGenMETPt");
  AddBranch(&patGenMETPhi_,    "patGenMETPhi");
  AddBranch(&patGenMETSumEt_,  "patGenMETSumEt");

  AddBranch(&patCaloMETPt_,     "patCaloMETPt");
  AddBranch(&patCaloMETPhi_,    "patCaloMETPhi");
  AddBranch(&patCaloMETSumEt_,  "patCaloMETSumEt");

  
  AddBranch(&puppiMETPt_,     "puppiMETPt");
  AddBranch(&puppiMETPhi_,    "puppiMETPhi");
  AddBranch(&puppiMETSumEt_,  "puppiMETSumEt");
  
  AddBranch(&puppiMETSig_,    "puppiMETSig");
  AddBranch(&puppiMETUnc_,    "puppiMETUnc");
  
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
  
  patGenMETPt_   = dummy ;
  patGenMETPhi_   = dummy ;
  patGenMETSumEt_   = dummy ;
  
  patCaloMETPt_   = dummy ;
  patCaloMETPhi_   = dummy ;
  patCaloMETSumEt_   = dummy ;

  puppiMETPt_     = dummy ;
  puppiMETPhi_    = dummy ;
  puppiMETSumEt_  = dummy ;
  
  puppiMETSig_    = dummy ;
  puppiMETUnc_.clear();
  
  
  


}
