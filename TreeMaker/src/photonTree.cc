
#include "ExoPieElement/TreeMaker/interface/photonTree.h"

photonTree::photonTree(std::string name, TTree* tree):
  baseTree(name,tree)
{
  photonP4_ =   new TClonesArray("TLorentzVector");
  SetBranches();
}

photonTree::~photonTree(){
  delete photonP4_;
}

void photonTree::Fill(const edm::Event& iEvent){
  Clear();
  //fetch the input collection
  //
  edm::Handle<edm::View<pat::Photon> > photonHandle;
  if(not iEvent.getByToken(photonToken,photonHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found: photon" <<std::endl; 
    exit(0);
  }  
  //pat::PhotonCollection phColl(*(photonHandle.product()));

 //sort the objects by transverse momentum
  //std::sort(phColl.begin(),phColl.end(),PtGreater());

  // IDs 
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;

  edm::Handle<edm::ValueMap<float> > mvaValues;
  
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;

  
  iEvent.getByToken(phoLooseIdMapToken,  loose_id_decisions);
  iEvent.getByToken(phoMediumIdMapToken,  medium_id_decisions);
  iEvent.getByToken(phoTightIdMapToken,  tight_id_decisions);
  iEvent.getByToken(phoMVAValuesMapToken, mvaValues);
  
  
  iEvent.getByToken(phoChargedIsolationToken,       phoChargedIsolationMap);
  iEvent.getByToken(phoNeutralHadronIsolationToken, phoNeutralHadronIsolationMap);
  iEvent.getByToken(phoPhotonIsolationToken,        phoPhotonIsolationMap);


  edm::View<pat::Photon>::const_iterator ph;
  //pat::PhotonCollection::const_iterator ph;
  for(ph=photonHandle->begin(); ph!=photonHandle->end(); ph++){
    if(ph->pt() < 14.) continue;
    if(TMath::Abs(ph->eta()) > 2.5) continue;
    nPho_++;
    new( (*photonP4_)[nPho_-1]) TLorentzVector(
					       ph->p4().px(),
					       ph->p4().py(),
					       ph->p4().pz(),
					       ph->p4().energy()
					       );
    //Ids
    const auto pho = photonHandle->ptrAt(nPho_-1);
    // std::cout<<" loose id = "<<(*loose_id_decisions)[pho]<<std::endl;

    isPassLoose.push_back((*loose_id_decisions)[pho]);
    isPassMedium.push_back((*medium_id_decisions)[pho]);
    isPassTight.push_back((*tight_id_decisions)[pho]);
    phoIDMVA_.push_back((*mvaValues)[pho]);


    phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho]);
    phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho]);
    phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho]);

    // --------
    // Photon ID variables used to make photon ID booleans. 
    // --------
    
    phoSCE_           .push_back((*ph).superCluster()->energy());
    phoSCRawE_        .push_back((*ph).superCluster()->rawEnergy());
    phoSCEta_         .push_back((*ph).superCluster()->eta());
    phoSCPhi_         .push_back((*ph).superCluster()->phi());
    phoSCEtaWidth_    .push_back((*ph).superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back((*ph).superCluster()->phiWidth());
    phoSCBrem_        .push_back((*ph).superCluster()->phiWidth()/(*ph).superCluster()->etaWidth());
    phohasPixelSeed_  .push_back((Int_t)ph->hasPixelSeed());
    phoEleVeto_       .push_back((Int_t)ph->passElectronVeto());
    phoR9_            .push_back(ph->r9());
    phoHoverE_        .push_back(ph->hadTowOverEm());
    phoSigmaIEtaIEta_ .push_back(ph->see());
    phoSigmaIEtaIPhi_ .push_back(ph->sep());
    phoSigmaIPhiIPhi_ .push_back(ph->spp());
    phoSigmaIEtaIEtaFull5x5_ .push_back(ph->full5x5_sigmaIetaIeta());
    phoR9Full5x5_            .push_back(ph->full5x5_r9());
    



  }  

}

void photonTree::SetBranches(){
  AddBranch(&nPho_  ,"nPho");
  AddBranch(&photonP4_, "phoP4");
  AddBranch(&isPassTight,"phoIsPassTight");
  AddBranch(&isPassLoose,"phoIsPassLoose");
  AddBranch(&isPassMedium,"phoIsPassMedium");
  AddBranch(&phoIDMVA_,"phoIDMVA");
  
  AddBranch(&phoSCE_,"phoSCE");
  AddBranch(&phoSCRawE_,"phoSCRawE");
  AddBranch(&phoSCEta_,"phoSCEta");
  AddBranch(&phoSCPhi_,"phoSCPhi");
  AddBranch(&phoSCEtaWidth_,"phoSCEtaWidth");
  AddBranch(&phoSCPhiWidth_,"phoSCPhiWidth");
  AddBranch(&phoSCBrem_,"phoSCBrem");
  AddBranch(&phohasPixelSeed_,"phohasPixelSeed");
  AddBranch(&phoEleVeto_,"phoEleVeto");
  AddBranch(&phoR9_,"phoR9");
  AddBranch(&phoHoverE_,"phoHoverE");
  AddBranch(&phoSigmaIEtaIEta_,"phoSigmaIEtaIEta");
  AddBranch(&phoSigmaIEtaIPhi_,"phoSigmaIEtaIPhi");
  AddBranch(&phoSigmaIPhiIPhi_,"phoSigmaIPhiIPhi");
  AddBranch(&phoSigmaIEtaIEtaFull5x5_,"phoSigmaIEtaIEtaFull5x5");
  AddBranch(&phoR9Full5x5_,"phoR9Full5x5");
  AddBranch(&phoPFChIso_,"phoPFChIso");
  AddBranch(&phoPFPhoIso_,"phoPFPhoIso");
  AddBranch(&phoPFNeuIso_,"phoPFNeuIso");

}

void photonTree::Clear(){
  nPho_ = 0; 
  photonP4_->Clear();
  isPassLoose.clear();
  isPassMedium.clear();
  isPassTight.clear();
  phoIDMVA_.clear();
  phoSCE_.clear();
  phoSCRawE_.clear();
  phoSCEta_.clear();
  phoSCPhi_.clear();
  phoSCEtaWidth_.clear();
  phoSCPhiWidth_.clear();
  phoSCBrem_.clear();
  phohasPixelSeed_.clear();
  
  phoEleVeto_.clear();
  phoR9_.clear();
  phoHoverE_.clear();
  phoSigmaIEtaIEta_.clear();
  phoSigmaIEtaIPhi_.clear();
  phoSigmaIPhiIPhi_.clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoR9Full5x5_.clear();

  phoPFChIso_.clear();
  phoPFPhoIso_.clear();
  phoPFNeuIso_.clear();
}
