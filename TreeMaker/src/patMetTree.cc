#include "ExoPieElement/TreeMaker/interface/patMetTree.h"
#include "ExoPieElement/TreeMaker/interface/eventInfo.h"
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
  is_Data = iEvent.isRealData();

  edm::Handle<pat::METCollection> patMetHandle;
  if(not iEvent.getByToken(pfMETToken,patMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"pfMet"<<std::endl; exit(0);}


  // adding modified Type-1 MET to the tree: EE Fixed
  edm::Handle<pat::METCollection> patMetModifiedHandle;
  if(not iEvent.getByToken(pfMETModifiedToken,patMetModifiedHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
	     <<"modified pfMet"<<std::endl; exit(0);}

  // adding PUPPI MET
  edm::Handle<pat::METCollection> puppiMetHandle;
  if(not iEvent.getByToken(puppimetToken, puppiMetHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: "
             <<"slimmedMETsPuppi"<<std::endl; exit(0);}




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
  if (!is_Data){
  patGenMETPt_ = met->genMET()->et();
  patGenMETPhi_ = met->genMET()->phi();
  patGenMETSumEt_ = met->genMET()->sumEt();
  }


  // calo met
  patCaloMETPt_ = met->caloMETPt();
  patCaloMETPhi_ = met->caloMETPhi();
  patCaloMETSumEt_ = met->caloMETSumEt();


  // CHS MET
  CHSMETPt_     = met->corPt(pat::MET::RawChs);
  CHSMETPhi_    = met->corPhi(pat::MET::RawChs);
  CHSMETSumEt_  = met->corSumEt(pat::MET::RawChs);

  // Track MET
  TRKMETPt_     = met->corPt(pat::MET::RawTrk);
  TRKMETPhi_    = met->corPhi(pat::MET::RawTrk);
  TRKMETPSumEt_ = met->corSumEt(pat::MET::RawTrk);

  // met uncertainties, need to be changed when time comes, right now don't have much time to edit this part
  patMet_smear_    = met->shiftedPt(pat::MET::NoShift, pat::MET::Type1Smear);

  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetResDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::JetEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::UnclusteredEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnUp));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::PhotonEnDown));
  patMetCorrUnc_.push_back(met->shiftedPt(pat::MET::NoShift));

  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::JetResUp));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::JetResDown));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::JetEnUp));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::JetEnDown));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::UnclusteredEnUp));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::UnclusteredEnDown));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::PhotonEnUp));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::PhotonEnDown));
  patMetCorrUnc_.push_back(met->shiftedPhi(pat::MET::NoShift));



  // Modified Type 1 corrected MET default in miniaod :: Needed only for 2017 data mc.

  pat::METCollection::const_iterator metmodified=patMetModifiedHandle.product()->begin();
  patmodifiedMetCorrPt_    = metmodified->et();
  patmodifiedMetCorrPhi_   = metmodified->phi();
  patmodifiedMetCorrSumEt_ = metmodified->sumEt();
  patmodifiedMetCorrSig_   = metmodified->significance() < 1.e10 ? met->significance() : 0;

  // met uncertainties, need to be changed when time comes, right now don't have much time to edit this part
  patmodifiedMet_smear_    = metmodified->shiftedPt(pat::MET::NoShift, pat::MET::Type1Smear);
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::JetResUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::JetResDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::JetEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::JetEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::UnclusteredEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::UnclusteredEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::PhotonEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::PhotonEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPt(pat::MET::NoShift));

  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::JetResUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::JetResDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::JetEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::JetEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::UnclusteredEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::UnclusteredEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::PhotonEnUp));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::PhotonEnDown));
  patmodifiedMetCorrUnc_.push_back(metmodified->shiftedPhi(pat::MET::NoShift));

// puppi met,  present in miniaod

  auto metpuppi = puppiMetHandle.product()->begin();
  puppiMETPt_         = metpuppi->et();
  puppiMETPhi_        = metpuppi->phi();
  puppiMETSumEt_      = metpuppi->sumEt();
  puppiMETSig_        = metpuppi->significance() < 1.e10 ? met->significance() : 0 ;


  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetResDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::JetEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::UnclusteredEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnUp));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::PhotonEnDown));
  puppiMETUnc_.push_back(metpuppi->shiftedPt(pat::MET::NoShift));

}
bool met_extra = false;
void
patMetTree::SetBranches(){

  AddBranch(&patMetCorrPt_, "MetCorrPt");
  AddBranch(&patMetCorrPhi_, "MetCorrPhi");
  AddBranch(&patMetCorrSumEt_, "MetCorrSumEt");
  AddBranch(&patMetCorrUnc_, "MetCorrUnc");
  AddBranch(&patMetCorrSig_, "MetCorrSig");

  AddBranch(&patMet_smear_, "patMet_smear");
  
  AddBranch(&patmodifiedMetCorrPt_, "modifiedMetCorrPt");
  AddBranch(&patmodifiedMetCorrPhi_, "modifiedMetCorrPhi");
  AddBranch(&patmodifiedMetCorrSumEt_, "modifiedMetCorrSumEt");
  AddBranch(&patmodifiedMetCorrSig_, "modifiedMetCorrSig");
  AddBranch(&patmodifiedMetCorrUnc_,"modifiedMetCorrUnc");

  AddBranch(&patmodifiedMet_smear_, "patmodifiedMet_smear");

  AddBranch(&patMetRawPt_, "MetRawPt");
  AddBranch(&patMetRawPhi_, "MetRawPhi");
  AddBranch(&patMetRawSumEt_, "MetRawSumEt");

  AddBranch(&patGenMETPt_,     "patGenMETPt");
  AddBranch(&patGenMETPhi_,    "patGenMETPhi");
  AddBranch(&patGenMETSumEt_,  "patGenMETSumEt");

  AddBranch(&patCaloMETPt_,     "patCaloMETPt");
  AddBranch(&patCaloMETPhi_,    "patCaloMETPhi");
  AddBranch(&patCaloMETSumEt_,  "patCaloMETSumEt");

  AddBranch(&CHSMETPt_, "CHSMETPt_");
  AddBranch(&CHSMETPhi_, "CHSMETPhi_");
  AddBranch(&CHSMETSumEt_, "CHSMETSumEt_");

  AddBranch(&TRKMETPt_, "TRKMETPt_");
  AddBranch(&TRKMETPhi_, "TRKMETPhi_");
  AddBranch(&TRKMETPSumEt_, "TRKMETPSumEt_");

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

  patMet_smear_= dummy;

  patmodifiedMetCorrPt_= dummy;
  patmodifiedMetCorrPhi_= dummy;
  patmodifiedMetCorrSumEt_= dummy;

  patmodifiedMetCorrSig_= dummy;
  patmodifiedMetCorrUnc_.clear();

  patmodifiedMet_smear_ =dummy;
  patMetRawPt_= dummy;
  patMetRawPhi_= dummy;
  patMetRawSumEt_= dummy;

  patGenMETPt_   = dummy ;
  patGenMETPhi_   = dummy ;
  patGenMETSumEt_   = dummy ;

  patCaloMETPt_   = dummy ;
  patCaloMETPhi_   = dummy ;
  patCaloMETSumEt_   = dummy ;

  CHSMETPt_   = dummy;
  CHSMETPhi_   = dummy;
  CHSMETSumEt_   = dummy;

  TRKMETPt_   = dummy;
  TRKMETPhi_   = dummy;
  TRKMETPSumEt_   = dummy;

  puppiMETPt_     = dummy ;
  puppiMETPhi_    = dummy ;
  puppiMETSumEt_  = dummy ;

  puppiMETSig_    = dummy ;
  puppiMETUnc_.clear();





}
