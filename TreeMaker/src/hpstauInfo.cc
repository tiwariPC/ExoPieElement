#include "ExoPieElement/TreeMaker/interface/hpstauInfo.h"

hpstauInfo::hpstauInfo(std::string name, TTree* tree, bool debug):baseTree(name,tree){
  if(debug) std::cout<<"in tau constructor"<<std::endl;
  HPSTau_4Momentum           = new TClonesArray("TLorentzVector");
  HPSTau_Vposition           = new TClonesArray("TVector3");
  if(debug) std::cout<<"in rho constructor: calling SetBrances()"<<std::endl;
  SetBranches();

  debug_ = debug ;
}

hpstauInfo::~hpstauInfo(){
  delete HPSTau_4Momentum;
  delete HPSTau_Vposition;
}

void hpstauInfo::Fill(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Clear();

  if(debug_)    std::cout<<"getting HPS Tau discriminators"<<std::endl;

    edm::Handle<pat::TauCollection>  tauHandle;
  if(not iEvent.getByToken(tauToken,tauHandle)){
    std::cout<<"FATAL EXCEPTION: "<<"Following Not Found:selectedPatTaus "<<std::endl;
    exit(0);
  }

  pat::TauCollection::const_iterator tau;
  for(tau=tauHandle->begin(); tau!=tauHandle->end(); tau++){
    if(tau->pt() < 18.0) continue;
    if(TMath::Abs(tau->eta()) > 2.5) continue;

    if(debug_) std::cout<<" pt,eta,phi of "<<HPSTau_n+1
			<<" tau is : "<<tau->pt()
			<<" : "<<tau->eta()
			<<" : "<<tau->phi()
			<<std::endl;
    TLorentzVector p4(tau->px(),tau->py(),tau->pz(),tau->energy());
    new( (*HPSTau_4Momentum)[HPSTau_n]) TLorentzVector(p4);
    TVector3 v3(tau->vx(),tau->vy(),tau->vz());
    new( (*HPSTau_Vposition)[HPSTau_n]) TVector3(v3);
    HPSTau_charge.push_back((int)tau->charge());

    TauPx_.push_back(tau->px());
    TauPy_.push_back(tau->py());
    TauPz_.push_back(tau->pz());
    TauE_.push_back(tau->energy());

    // New Tau Disc for Fall17 V2 
    
    // -- decay mode finding 
    disc_decayModeFinding.push_back(tau->tauID("decayModeFinding"));
    disc_decayModeFindingNewDMs.push_back(tau->tauID("decayModeFindingNewDMs"));
    
    // -- isolation 
    disc_byIsolationMVArun2017v2DBoldDMwLTraw2017.push_back(tau->tauID("byIsolationMVArun2017v2DBoldDMwLTraw2017"));
    disc_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byVLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byLooseIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byMediumIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byTightIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byVTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVTightIsolationMVArun2017v2DBoldDMwLT2017"));
    disc_byVVTightIsolationMVArun2017v2DBoldDMwLT2017.push_back(tau->tauID("byVVTightIsolationMVArun2017v2DBoldDMwLT2017"));
    
    // -- against electron 
    
    disc_againstElectronLooseMVA6.push_back(tau->tauID("againstElectronLooseMVA6"));
    disc_againstElectronMediumMVA6.push_back(tau->tauID("againstElectronMediumMVA6"));
    disc_againstElectronTightMVA6.push_back(tau->tauID("againstElectronTightMVA6"));
    disc_againstElectronVLooseMVA6.push_back(tau->tauID("againstElectronVLooseMVA6"));
    disc_againstElectronVTightMVA6.push_back(tau->tauID("againstElectronVTightMVA6"));
    
    // -- against muons

    disc_againstMuonLoose3.push_back(tau->tauID("againstMuonLoose3"));
    disc_againstMuonTight3.push_back(tau->tauID("againstMuonTight3"));

    // -- cut based isolation   
    disc_byLooseCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
    disc_byMediumCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
    disc_byTightCombinedIsolationDeltaBetaCorr3Hits.push_back(tau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

    // -- pt sum for optimisation, going to switch off them, 
    disc_chargedIsoPtSum.push_back(tau->tauID("chargedIsoPtSum"));
    disc_neutralIsoPtSum.push_back(tau->tauID("neutralIsoPtSum"));
    disc_puCorrPtSum.push_back(tau->tauID("puCorrPtSum"));


    HPSTau_leadPFChargedHadrCand.push_back(tau->leadPFChargedHadrCand().isNonnull() );
    HPSTau_leadPFChargedHadrCand_trackRef.push_back(tau->leadPFChargedHadrCand().isNonnull() &&  tau->leadPFChargedHadrCand()->trackRef().isNonnull());

    //std::cout<<" middle pt = "<<tau->pt()<<std::endl;

    using namespace reco;
    using namespace pat;
    using namespace edm;

    ESHandle<MagneticField> B;
    iSetup.get<IdealMagneticFieldRecord > ().get(B);
    const MagneticField* magField = B.product();
    Handle<reco::BeamSpot> theBeamSpotHandle;
    iEvent.getByToken(theBeamSpotToken, theBeamSpotHandle);
    const reco::BeamSpot* beamSpot = theBeamSpotHandle.product();
    ESHandle<GlobalTrackingGeometry> geomHandle;
    iSetup.get<GlobalTrackingGeometryRecord > ().get(geomHandle);

    float newvz=0.0;
    if(tau->leadPFChargedHadrCand().isNonnull()) {
      if (tau->leadPFChargedHadrCand()->trackRef().isNonnull()) {
	reco::TransientTrack track(tau->leadPFChargedHadrCand()->trackRef(), magField, geomHandle);
	TransverseImpactPointExtrapolator extrapolator(magField);
	TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(), GlobalPoint(beamSpot->position().x(), beamSpot->position().y(), 0.0));
	newvz = (closestOnTransversePlaneState.globalPosition().z());
      }

      if(tau->leadPFChargedHadrCand()->gsfTrackRef().isNonnull()) {
	reco::GsfTransientTrack track(tau->leadPFChargedHadrCand()->gsfTrackRef(),magField,geomHandle);
	TransverseImpactPointExtrapolator extrapolator(magField);
	TrajectoryStateOnSurface closestOnTransversePlaneState = extrapolator.extrapolate(track.impactPointState(),GlobalPoint(beamSpot->position().x(),beamSpot->position().y(),0.0));
	newvz = (closestOnTransversePlaneState.globalPosition().z());
      }
    }

    HPSTau_NewVz.push_back(newvz);
    taupt.push_back(tau->pt());
    //std::cout<<"other  pt = "<<tau->pt()<<"  :  new vz = "<<newvz<<std::endl;
    HPSTau_n++;
  }//end of for loop
  // std::cout<<" -------------------------- HPSTau_n = "<<HPSTau_n<<std::endl;
  if(debug_)    std::cout<<"got HPS Tau discriminators info"<<std::endl;
}



void hpstauInfo::SetBranches(){
  if(debug_)    std::cout<<"setting branches: calling AddBranch of baseTree"<<std::endl;

  AddBranch(&HPSTau_n  ,"HPSTau_n");
  AddBranch(&taupt  ,"taupt");

//  AddBranch(&HPSTau_4Momentum,"HPSTau_4Momentum");
  AddBranch(&HPSTau_Vposition,"HPSTau_Vposition");

  AddBranch(&TauPx_, "HPSTau_Px");
  AddBranch(&TauPy_, "HPSTau_Py");
  AddBranch(&TauPz_, "HPSTau_Pz");
  AddBranch(&TauE_, "HPSTau_Energy");

  if (false){ 
    // these are all old, just keeping them for now ,will be killed in future, 
    AddBranch(&HPSTau_leadPFChargedHadrCand,"HPSTau_leadPFChargedHadrCand");
    AddBranch(&HPSTau_leadPFChargedHadrCand_trackRef,"HPSTau_leadPFChargedHadrCand_trackRef");
    
    
    AddBranch(&disc_againstElectronLoose ,"disc_againstElectronLoose");
    AddBranch(&disc_againstElectronMedium ,"disc_againstElectronMedium");
    AddBranch(&disc_againstElectronTight ,"disc_againstElectronTight");
    
    AddBranch(&disc_againstMuonLoose ,"disc_againstMuonLoose");
    AddBranch(&disc_againstMuonMedium ,"disc_againstMuonMedium");
    AddBranch(&disc_againstMuonTight ,"disc_againstMuonTight");
    AddBranch(&disc_againstMuonLoose2 ,"disc_againstMuonLoose2");
    AddBranch(&disc_againstMuonMedium2 ,"disc_againstMuonMedium2");
    AddBranch(&disc_againstMuonTight2 ,"disc_againstMuonTight2");
    AddBranch(&disc_againstMuonLooseMVA ,"disc_againstMuonLooseMVA");
    AddBranch(&disc_againstMuonMediumMVA ,"disc_againstMuonMediumMVA");
    AddBranch(&disc_againstMuonTightMVA ,"disc_againstMuonTightMVA");
    AddBranch(&disc_byVLooseCombinedIsolationDeltaBetaCorr ,"disc_byVLooseCombinedIsolationDeltaBetaCorr");
    AddBranch(&disc_byLooseIsolation ,"disc_byLooseIsolation");
    // 2017 v2 new DM
    AddBranch(&disc_byIsolationMVArun2017v2DBnewDMwLTraw2017,"disc_byIsolationMVArun2017v2DBnewDMwLTraw2017");
    AddBranch(&disc_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017,"disc_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byVLooseIsolationMVArun2017v2DBnewDMwLT2017,"disc_byVLooseIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byLooseIsolationMVArun2017v2DBnewDMwLT2017,"disc_byLooseIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byMediumIsolationMVArun2017v2DBnewDMwLT2017,"disc_byMediumIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byTightIsolationMVArun2017v2DBnewDMwLT2017,"disc_byTightIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byVTightIsolationMVArun2017v2DBnewDMwLT2017,"disc_byVTightIsolationMVArun2017v2DBnewDMwLT2017");
    AddBranch(&disc_byVVTightIsolationMVArun2017v2DBnewDMwLT2017,"disc_byVVTightIsolationMVArun2017v2DBnewDMwLT2017");

    AddBranch(&disc_chargedIsoPtSum ,"disc_chargedIsoPtSum");
    AddBranch(&disc_neutralIsoPtSum ,"disc_neutralIsoPtSum");
    AddBranch(&disc_puCorrPtSum ,"disc_puCorrPtSum");

  }
  
  AddBranch(&disc_againstElectronLooseMVA6 ,"disc_againstElectronLooseMVA6");  
  AddBranch(&disc_againstElectronMediumMVA6 ,"disc_againstElectronMediumMVA6");
  AddBranch(&disc_againstElectronTightMVA6 ,"disc_againstElectronTightMVA6");
  AddBranch(&disc_againstElectronVLooseMVA6 ,"disc_againstElectronVLooseMVA6");
  AddBranch(&disc_againstElectronVTightMVA6 ,"disc_againstElectronVTightMVA6");
    
  
  AddBranch(&disc_againstMuonLoose3 ,"disc_againstMuonLoose3");
  AddBranch(&disc_againstMuonTight3 ,"disc_againstMuonTight3");

  AddBranch(&disc_byLooseCombinedIsolationDeltaBetaCorr ,"disc_byLooseCombinedIsolationDeltaBetaCorr");
  AddBranch(&disc_byMediumCombinedIsolationDeltaBetaCorr ,"disc_byMediumCombinedIsolationDeltaBetaCorr");
  AddBranch(&disc_byTightCombinedIsolationDeltaBetaCorr ,"disc_byTightCombinedIsolationDeltaBetaCorr");

  

  AddBranch(&disc_byIsolationMVArun2017v2DBoldDMwLTraw2017,"disc_byIsolationMVArun2017v2DBoldDMwLTraw2017");
  AddBranch(&disc_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017,"disc_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&disc_byVLooseIsolationMVArun2017v2DBoldDMwLT2017,"disc_byVLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017,"disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017,"disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(&disc_byTightIsolationMVArun2017v2DBoldDMwLT2017,"disc_byTightIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(& disc_byVTightIsolationMVArun2017v2DBoldDMwLT2017,"disc_byVTightIsolationMVArun2017v2DBoldDMwLT2017");
  AddBranch(& disc_byVVTightIsolationMVArun2017v2DBoldDMwLT2017,"disc_byVVTightIsolationMVArun2017v2DBoldDMwLT2017");
  
  

  AddBranch(&disc_byLooseCombinedIsolationDeltaBetaCorr3Hits ,"disc_byLooseCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&disc_byMediumCombinedIsolationDeltaBetaCorr3Hits ,"disc_byMediumCombinedIsolationDeltaBetaCorr3Hits");
  AddBranch(&disc_byTightCombinedIsolationDeltaBetaCorr3Hits ,"disc_byTightCombinedIsolationDeltaBetaCorr3Hits");

  AddBranch(&disc_decayModeFinding ,"disc_decayModeFinding");
  AddBranch(&disc_decayModeFindingNewDMs ,"disc_decayModeFindingNewDMs");



  AddBranch(&HPSTau_NewVz,"HPSTau_NewVz");
  AddBranch(&HPSTau_charge,"HPSTau_charge");

  if(debug_)    std::cout<<"set branches"<<std::endl;
}

void hpstauInfo::Clear(){
  if(debug_)std::cout<<"clearing HpsTau info"<<std::endl;
  HPSTau_n =0;
  HPSTau_4Momentum->Clear();
  HPSTau_Vposition->Clear();
  taupt.clear();

  TauPx_.clear();
  TauPy_.clear();
  TauPz_.clear();
  TauE_.clear();

  HPSTau_leadPFChargedHadrCand.clear();
  HPSTau_leadPFChargedHadrCand_trackRef.clear();



  disc_againstElectronLoose.clear();
  disc_againstElectronMedium.clear();
  disc_againstElectronTight.clear();
  disc_againstElectronLooseMVA6.clear();
  disc_againstElectronMediumMVA6.clear();
  disc_againstElectronTightMVA6.clear();
  disc_againstElectronVLooseMVA6.clear();
  disc_againstElectronVTightMVA6.clear();
  disc_againstMuonLoose.clear();
  disc_againstMuonMedium.clear();
  disc_againstMuonTight.clear();
  disc_againstMuonLoose2.clear();
  disc_againstMuonMedium2.clear();
  disc_againstMuonTight2.clear();
  disc_againstMuonLooseMVA.clear();
  disc_againstMuonMediumMVA.clear();
  disc_againstMuonTightMVA.clear();
  disc_againstMuonLoose3.clear();
  disc_againstMuonTight3.clear();

  disc_byVLooseCombinedIsolationDeltaBetaCorr.clear();
  disc_byLooseCombinedIsolationDeltaBetaCorr.clear();
  disc_byMediumCombinedIsolationDeltaBetaCorr.clear();
  disc_byTightCombinedIsolationDeltaBetaCorr.clear();

  disc_byLooseIsolation.clear();


  disc_byIsolationMVArun2017v2DBoldDMwLTraw2017.clear();
  disc_byVVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byVLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byLooseIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byMediumIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byVTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  disc_byVVTightIsolationMVArun2017v2DBoldDMwLT2017.clear();
  // 2017 v2 new DM
  disc_byIsolationMVArun2017v2DBnewDMwLTraw2017.clear();
  disc_byVVLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byVLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byLooseIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byMediumIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byTightIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byVTightIsolationMVArun2017v2DBnewDMwLT2017.clear();
  disc_byVVTightIsolationMVArun2017v2DBnewDMwLT2017.clear();

  disc_byLooseCombinedIsolationDeltaBetaCorr3Hits.clear();
  disc_byMediumCombinedIsolationDeltaBetaCorr3Hits.clear();
  disc_byTightCombinedIsolationDeltaBetaCorr3Hits.clear();

  disc_decayModeFinding.clear();
  disc_decayModeFindingNewDMs.clear();

  disc_chargedIsoPtSum.clear();
  disc_neutralIsoPtSum.clear();
  disc_puCorrPtSum.clear();

  HPSTau_NewVz.clear();
  HPSTau_charge.clear();

  if(debug_) std::cout<<"cleared"<<std::endl;
}
