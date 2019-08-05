#include "ExoPieElement/TreeMaker/interface/patElecTree.h"
// for conversion finder
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

// Issues to be resolved :
// -- mypv
// -- rho
// -- fix the impact parameter
// -- make conv veto boolean
// -- fix the initial values for isolation

patElecTree::patElecTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(name,tree),
  r_iso_min_(iConfig.getParameter<double>("r_iso_min")),
  r_iso_max_(iConfig.getParameter<double>("r_iso_max")),
  kt_scale_(iConfig.getParameter<double>("kt_scale")),
  charged_only_(iConfig.getParameter<bool>("charged_only")),
  eAreasElectrons("effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt")
{
  patElecP4_ =   new TClonesArray("TLorentzVector");
  SetBranches();

}
patElecTree::~patElecTree(){
  delete patElecP4_;
}

void
patElecTree::Fill(const edm::Event& iEvent){
  Clear();


  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByToken(eleToken,electronHandle);

  /*  //id boolean
  

  edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
  edm::Handle<edm::ValueMap<bool> > heep_id_decisions;

  edm::Handle<ValueMap<vid::CutFlowResult>> veto_id_cutflow;
  edm::Handle<ValueMap<vid::CutFlowResult>> loose_id_cutflow;
  edm::Handle<ValueMap<vid::CutFlowResult>> medium_id_cutflow;
  edm::Handle<ValueMap<vid::CutFlowResult>> tight_id_cutflow;
  edm::Handle<ValueMap<vid::CutFlowResult>> heep_id_cutflow;

  edm::Handle<edm::ValueMap<bool> > medium_MVAid_decisions;
  edm::Handle<edm::ValueMap<bool> > tight_MVAid_decisions;

  std::cout<< "just before ele id decision "<<std::endl;
  
  iEvent.getByToken(eleVetoIdMapToken,   veto_id_decisions);
  iEvent.getByToken(eleVetoIdCFToken,    veto_id_cutflow);


  std::cout<<" just before loose "<<std::endl;

  iEvent.getByToken(eleLooseIdMapToken,  loose_id_decisions);
  iEvent.getByToken(eleLooseIdCFToken,   loose_id_cutflow);

  std::cout<<" just before medium"<<std::endl;
    
  iEvent.getByToken(eleMediumIdMapToken, medium_id_decisions);
  iEvent.getByToken(eleMediumIdCFToken,  medium_id_cutflow);

  iEvent.getByToken(eleTightIdMapToken,  tight_id_decisions);
  iEvent.getByToken(eleTightIdCFToken,   tight_id_cutflow);

  std::cout<<" just before heep"<<std::endl;

  iEvent.getByToken(eleHEEPIdMapToken,   heep_id_decisions);
  iEvent.getByToken(eleHEEPIdCFToken,    heep_id_cutflow);

  std::vector<std::string> maskCutBasedCuts;
  maskCutBasedCuts.push_back("GsfEleEffAreaPFIsoCut_0");
  std::vector<std::string> maskHEEPCuts;
  //maskHEEPCuts.push_back("GsfEleTrkPtIsoCut_0");
  maskHEEPCuts.push_back("GsfEleEmHadD1IsoRhoCut_0");

  iEvent.getByToken(eleMVAMediumIdMapToken, medium_MVAid_decisions);
  iEvent.getByToken(eleMVATightIdMapToken,  tight_MVAid_decisions);

  // Get MVA values and categories (optional)
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<int> > mvaCategories;
  iEvent.getByToken(mvaValuesMapToken,      mvaValues);
  iEvent.getByToken(mvaCategoriesMapToken,  mvaCategories);

  */

  edm::Handle<reco::VertexCollection> recVtxs;
  if(not iEvent.getByToken(vertexToken, recVtxs))return;

  if (recVtxs->empty()) return; // skip the event if no PV found
  vector<reco::Vertex>::const_iterator firstGoodVertex = recVtxs->end();
  //VertexCollection::const_iterator firstGoodVertex = recVtxs->end();

  std::cout<<" just before vertex loop"<<std::endl;
  int firstGoodVertexIdx = 0;
  //  for (VertexCollection::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx, ++firstGoodVertexIdx) {
  for (vector<reco::Vertex>::const_iterator vtx = recVtxs->begin(); vtx != recVtxs->end(); ++vtx, ++firstGoodVertexIdx) {

    // Replace isFake() for miniAOD because it requires tracks and miniAOD recVtxs don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    // bool isFake = vtx->isFake();
    //if( !isAOD ) //we are here for MINIAOD only
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);
    // Check the goodness
    if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      firstGoodVertex = vtx;
      break;
    }
  }

  if ( firstGoodVertex==recVtxs->end() )
    return; // skip event if there are no good PVs


  // handle pfcandidates
  Handle<pat::PackedCandidateCollection> pfcands;
  iEvent.getByToken(pfCandToken, pfcands);

  // Get rho value
  edm::Handle<double> rhoH;
  iEvent.getByToken(rhoForLepToken,rhoH);
  patElecRho_ = *rhoH;

  for (edm::View<pat::Electron>::const_iterator ele = electronHandle->begin(); ele != electronHandle->end(); ++ele) {

    if(ele->pt() < 10.) continue;
    if(TMath::Abs(ele->eta()) > 2.5) continue;
    nEle_++;
    
    std::cout<<" just before ele p4"<<std::endl;
    new( (*patElecP4_)[nEle_-1]) TLorentzVector(
						ele->p4().px(),
						ele->p4().py(),
						ele->p4().pz(),
						ele->p4().energy()
						);
//   px,py,pz,e
    patElecPx_.push_back(ele->p4().px());
    patElecPy_.push_back(ele->p4().py());
    patElecPz_.push_back(ele->p4().pz());
    patElecE_.push_back(ele->p4().energy());

    patElecInBarrel_.push_back(ele->isEB());
    patElecInEndcap_.push_back(ele->isEE());

    patElecCharge_.push_back(ele->charge());
    patElecChargeConsistent_.push_back(ele->isGsfCtfScPixChargeConsistent());

    std::cout<<" just before caloenergy"<<std::endl;
    patElecaloEnergy_.push_back(ele->caloEnergy());

    double R = sqrt(ele->superCluster()->x()*ele->superCluster()->x() + ele->superCluster()->y()*ele->superCluster()->y() +ele->superCluster()->z()*ele->superCluster()->z());
    double Rt = sqrt(ele->superCluster()->x()*ele->superCluster()->x() + ele->superCluster()->y()*ele->superCluster()->y());
    patElecScEt_.push_back( (ele->superCluster()->energy())*(Rt/R) );
    patElecScEn_.push_back(ele->superCluster()->energy());
    patElecScPreEn_.push_back(ele->superCluster()->preshowerEnergy());
    patElecScEta_.push_back(ele->superCluster()->eta());
    patElecScPhi_.push_back(ele->superCluster()->phi());
    patElecScRawEn_.push_back(ele->superCluster()->rawEnergy());
    patElecScEtaWidth_.push_back(ele->superCluster()->etaWidth());
    patElecScPhiWidth_.push_back(ele->superCluster()->phiWidth());

    patElecR9_.push_back(ele->r9());
    patElecHoverE_.push_back(ele->hcalOverEcal());

    // fix this
    patElecD0_.push_back(ele->gsfTrack()->dxy(firstGoodVertex->position()));
    patElecDz_.push_back(ele->gsfTrack()->dz(firstGoodVertex->position()));

    patElecEoverP_.push_back(ele->eSuperClusterOverP());
    patElecBrem_.push_back(ele->fbrem());
    patElecdEtaAtVtx_.push_back(ele->deltaEtaSuperClusterTrackAtVtx());
    patElecdPhiAtVtx_.push_back(ele->deltaPhiSuperClusterTrackAtVtx());
    patElecSigmaIEtaIEta_.push_back(ele->sigmaIetaIeta()); ///new sigmaietaieta
    patElecSigmaIEtaIPhi_.push_back(ele->sigmaIetaIphi());
    patElecSigmaIPhiIPhi_.push_back(ele->sigmaIphiIphi());


    std::cout<<" just before ele conv veto"<<std::endl;
    
    patElecConvVeto_.push_back(ele->passConversionVeto()); // ConvVtxFit || missHit == 0
    patElecMissHits_.push_back(ele->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
    if (ele->ecalEnergy() == 0) {
      patElecEoverPInv_.push_back(1e30);
    } else if (!std::isfinite(ele->ecalEnergy())) {
      patElecEoverPInv_.push_back(1e30);
    } else {
      patElecEoverPInv_.push_back(fabs(1./ele->ecalEnergy() - ele->eSuperClusterOverP()/ele->ecalEnergy()));
    }
    ///HEEP ID
    // double eledEtaseedAtVtx = ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    //   ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
    // patElecdEtaseedAtVtx_.push_back(eledEtaseedAtVtx);
    patElecdEtaseedAtVtx_.push_back(ele->deltaEtaSeedClusterTrackAtVtx());
    patElecE1x5_.push_back(ele->e1x5());
    patElecE2x5_.push_back(ele->e2x5Max());
    patElecE5x5_.push_back(ele->e5x5());

    // this variable is for debugging

    patElecSigmaIEtaIEtaFull5x5_.push_back(ele->full5x5_sigmaIetaIeta());
    patElecE1x5Full5x5_.push_back(ele->full5x5_e1x5());
    patElecE2x5Full5x5_.push_back(ele->full5x5_e2x5Max());
    patElecE5x5Full5x5_.push_back(ele->full5x5_e5x5());
    patElecR9Full5x5_.push_back(ele->full5x5_r9());
    

    std::cout<<" just before isolation vars"<<std::endl;
    //To include in anlyzer code
    /*    edm::FileInPath eaConstantsFile("EgammaAnalysis/ElectronTools/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt");
	  EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());
	  float eA = effectiveAreas.getEffectiveArea(fabs(ele->superCluster()->eta()));
	  patElecEffArea_.push_back(eA);    */

    //fix the initial value
    float iso1 = 999.;
    float iso2 = 999.;
    float iso3 = 999.;
    float isoPU = -999.;
    iso1 = ele->pfIsolationVariables().sumChargedHadronPt;
    iso2 = ele->pfIsolationVariables().sumNeutralHadronEt;
    iso3 = ele->pfIsolationVariables().sumPhotonEt;
    isoPU = ele->pfIsolationVariables().sumPUPt;
    patElecChHadIso_.push_back(iso1);
    patElecNeHadIso_.push_back(iso2);
    patElecGamIso_.push_back(iso3);
    patElecPUPt_.push_back(isoPU);

    
    /*
    double miniIso[7]={0};
    getPFIsolation(miniIso, pfcands, dynamic_cast<const reco::Candidate *>(&(*ele)),
		   eAreasElectrons, ele->superCluster()->eta(),
		   *rhoH, r_iso_min_, r_iso_max_, kt_scale_, charged_only_);

    std::cout<<" just before miniIso"<<std::endl;
    patElecMiniIso_ch_.push_back(miniIso[0]);
    patElecMiniIso_nh_.push_back(miniIso[1]);
    patElecMiniIso_ph_.push_back(miniIso[2]);
    patElecMiniIso_pu_.push_back(miniIso[3]);
    patElecMiniIso_r_.push_back(miniIso[4]);
    patElecMiniIsoBeta_.push_back(miniIso[5]);
    patElecMiniIsoEA_.push_back(miniIso[6]);

    ///For HEEP ID
    patElecEcalDrivenSeed_.push_back(ele->ecalDrivenSeed());
    patElecEcalDriven_.push_back(ele->ecalDriven());
    patElecDr03EcalRecHitSumEt_.push_back(ele->dr03EcalRecHitSumEt());
    patElecDr03HcalDepth1TowerSumEt_.push_back(ele->dr03HcalDepth1TowerSumEt());
    patElecDr03HcalDepth2TowerSumEt_.push_back(ele->dr03HcalDepth2TowerSumEt());
    patElecDr03HcalTowerSumEt_ .push_back(ele->dr03HcalTowerSumEt());
    patElecDr03TkSumPt_ .push_back(ele->dr03TkSumPt());


    // for MVA preselection
    patElecEcalPFClusterIso_.push_back(ele->ecalPFClusterIso());
    patElecHcalPFClusterIso_.push_back(ele->hcalPFClusterIso());

    */
    // reco::GsfTrackRef trackref = ele->gsfTrack();
    // Fix this impact parameter
    /* if (ele->gsfTrack().isNonnull()) {
       if (recVtxs->size() > 0){
       patElecTrkdz_.push_back(trackref->dz(recVtxs->front().position()));
       patElecTrkdxy_.push_back(trackref->dxy(recVtxs->front().position()));
       }
       }*/

    std::cout<<" just before saving id decision"<<std::endl;
    const auto el = electronHandle->ptrAt(nEle_-1);


    std::cout<<" just before saving veto id decision"<<std::endl;
    isPassVeto_.push_back(ele->electronID("cutBasedElectronID-Fall17-94X-V2-veto"));
    isPassLoose_.push_back(ele->electronID("cutBasedElectronID-Fall17-94X-V2-loose"));
    isPassMedium_.push_back(ele->electronID("cutBasedElectronID-Fall17-94X-V2-medium"));
    isPassTight_.push_back(ele->electronID("cutBasedElectronID-Fall17-94X-V2-tight"));
    isPassHEEP_.push_back(ele->electronID("heepElectronID-HEEPV70"));

    
  }
}
bool ele_extra = false;

void
patElecTree::SetBranches(){

  AddBranch(&patElecRho_, "eleRho");
  AddBranch(&nEle_, "nEle");


  AddBranch(&patElecPx_, "elePx");
  AddBranch(&patElecPy_, "elePy");
  AddBranch(&patElecPz_, "elePz");
  AddBranch(&patElecE_, "eleEnergy");

  AddBranch(&patElecCharge_, "eleCharge");
  AddBranch(&isPassVeto_,"eleIsPassVeto");
  AddBranch(&isPassLoose_,"eleIsPassLoose");
  AddBranch(&isPassMedium_,"eleIsPassMedium");
  AddBranch(&isPassTight_,"eleIsPassTight");
  
  //AddBranch(&isPassHEEP_,"eleIsPassHEEP");

  if (ele_extra){
    AddBranch(&patElecP4_,"eleP4");
    AddBranch(&patElecChargeConsistent_,"eleChargeConsistent");
    AddBranch(&patElecInBarrel_,"eleInBarrel");
    AddBranch(&patElecInEndcap_,"eleInEndcap");

    AddBranch(&patElecaloEnergy_,"elecaloEnergy");

    AddBranch(&patElecScEt_,"eleScEt");
    AddBranch(&patElecScEn_,"eleScEn");
    AddBranch(&patElecScPreEn_,"eleScPreEn");
    AddBranch(&patElecScEta_,"eleScEta");
    AddBranch(&patElecScPhi_,"eleScPhi");
    AddBranch(&patElecScRawEn_,"eleScRawEn");
    AddBranch(&patElecScEtaWidth_,"eleScEtaWidth");
    AddBranch(&patElecScPhiWidth_,"eleScPhiWidth");

    AddBranch(&patElecR9_,"eleR9");
    AddBranch(&patElecHoverE_,"eleHoverE");

    AddBranch(&patElecD0_,"eleD0");
    AddBranch(&patElecDz_,"eleDz");

    AddBranch(&patElecEoverP_,"eleEoverP");
    AddBranch(&patElecBrem_,"eleBrem");
    AddBranch(&patElecdEtaAtVtx_,"eledEtaAtVtx");
    AddBranch(&patElecdPhiAtVtx_,"eledPhiAtVtx");
    AddBranch(&patElecSigmaIEtaIEta_,"eleSigmaIEtaIEta");
    AddBranch(&patElecSigmaIEtaIPhi_,"eleSigmaIEtaIPhi");
    AddBranch(&patElecSigmaIPhiIPhi_,"eleSigmaIPhiIPhi");

    AddBranch(&patElecConvVeto_,"eleConvVeto");
    AddBranch(&patElecMissHits_,"eleMissHits");
    AddBranch(&patElecEoverPInv_,"eleEoverPInv");

    AddBranch(&patElecdEtaseedAtVtx_,"eleEtaseedAtVtx");
    AddBranch(&patElecE1x5_,"eleE1x5");
    AddBranch(&patElecE2x5_,"eleE2x5");
    AddBranch(&patElecE5x5_,"eleE5x5");

    AddBranch(&patElecSigmaIEtaIEtaFull5x5_,"eleSigmaIEtaIEtaFull5x5");
    AddBranch(&patElecE1x5Full5x5_,"eleE1x5Full5x5");
    AddBranch(&patElecE2x5Full5x5_,"eleE2x5Full5x5");
    AddBranch(&patElecE5x5Full5x5_,"eleE5x5Full5x5");
    AddBranch(&patElecR9Full5x5_,"eleR9Full5x5");

    AddBranch(&patElecChHadIso_, "eleChHadIso");
    AddBranch(&patElecNeHadIso_, "eleNeHadIso");
    AddBranch(&patElecGamIso_, "eleGamIso");
    AddBranch(&patElecPUPt_, "elePUPt");
    AddBranch(&patElecEcalPFClusterIso_, "eleEcalPFClusterIso");
    AddBranch(&patElecHcalPFClusterIso_, "eleHcalPFClusterIso");

    AddBranch(&patElecMiniIso_ch_,"eleMiniIso_ch");
    AddBranch(&patElecMiniIso_nh_,"eleMiniIso_nh");
    AddBranch(&patElecMiniIso_ph_,"eleMiniIso_ph");
    AddBranch(&patElecMiniIso_pu_,"eleMiniIso_pu");
    AddBranch(&patElecMiniIso_r_,"eleMiniIso_r");
    AddBranch(&patElecMiniIsoBeta_,"eleMiniIsoBeta");
    AddBranch(&patElecMiniIsoEA_,"eleMiniIsoEA");

    AddBranch(&patElecEcalDrivenSeed_,"eleEcalDrivenSeed");
    AddBranch(&patElecEcalDriven_,"eleEcalDriven");
    AddBranch(&patElecDr03EcalRecHitSumEt_,"eleDr03EcalRecHitSumEt");
    AddBranch(&patElecDr03HcalDepth1TowerSumEt_,"eleDr03HcalDepth1TowerSumEt");
    AddBranch(&patElecDr03HcalDepth2TowerSumEt_,"eleDr03HcalDepth2TowerSumEt");
    AddBranch(&patElecDr03HcalTowerSumEt_,"eleDr03HcalTowerSumEt");
    AddBranch(&patElecDr03TkSumPt_,"eleDr03TkSumPt");


    AddBranch(&isPassVetoNoIso_,"eleIsPassVetoNoIso");
    AddBranch(&isPassLooseNoIso_,"eleIsPassLooseNoIso");
    AddBranch(&isPassMediumNoIso_,"eleIsPassMediumNoIso");
    AddBranch(&isPassTightNoIso_,"eleIsPassTightNoIso");
    AddBranch(&isPassHEEPNoIso_,"eleIsPassHEEPNoIso");
    AddBranch(&isPassMVAMedium_,"eleIsPassMVAMedium");
    AddBranch(&isPassMVATight_,"eleIsPassMVATight");

    AddBranch(&mvaValue_,"eleMVAValue");
    AddBranch(&mvaCategory_,"eleMVACategory");
  }
}
void
patElecTree::Clear(){


  patElecRho_ =-999;
  nEle_ =0;
  patElecP4_->Clear();

  patElecPx_.clear();
  patElecPy_.clear();
  patElecPz_.clear();
  patElecE_.clear();

  patElecInBarrel_.clear();
  patElecInEndcap_.clear();


  patElecCharge_.clear();
  patElecChargeConsistent_.clear();

  patElecaloEnergy_.clear();

  patElecScEt_.clear();
  patElecScEn_.clear();
  patElecScPreEn_.clear();
  patElecScEta_.clear();
  patElecScPhi_.clear();
  patElecScRawEn_.clear();
  patElecScEtaWidth_.clear();
  patElecScPhiWidth_.clear();

  patElecR9_.clear();
  patElecHoverE_.clear();

  patElecD0_.clear();
  patElecDz_.clear();

  patElecEoverP_.clear();
  patElecBrem_.clear();
  patElecdEtaAtVtx_.clear();
  patElecdPhiAtVtx_.clear();
  patElecSigmaIEtaIEta_.clear();
  patElecSigmaIEtaIPhi_.clear();
  patElecSigmaIPhiIPhi_.clear();

  patElecConvVeto_.clear();
  patElecMissHits_.clear();
  patElecEoverPInv_.clear();

  patElecdEtaseedAtVtx_.clear();
  patElecE1x5_.clear();
  patElecE2x5_.clear();
  patElecE5x5_.clear();

  patElecSigmaIEtaIEtaFull5x5_.clear();
  patElecE1x5Full5x5_.clear();
  patElecE2x5Full5x5_.clear();
  patElecE5x5Full5x5_.clear();
  patElecR9Full5x5_.clear();

  patElecChHadIso_.clear();
  patElecNeHadIso_.clear();
  patElecGamIso_.clear();
  patElecPUPt_.clear();
  patElecEcalPFClusterIso_.clear();
  patElecHcalPFClusterIso_.clear();
  patElecMiniIso_ch_.clear();
  patElecMiniIso_nh_.clear();
  patElecMiniIso_ph_.clear();
  patElecMiniIso_pu_.clear();
  patElecMiniIso_r_.clear();
  patElecMiniIsoBeta_.clear();
  patElecMiniIsoEA_.clear();

  patElecEcalDrivenSeed_.clear();
  patElecEcalDriven_.clear();
  patElecDr03EcalRecHitSumEt_.clear();
  patElecDr03HcalDepth1TowerSumEt_.clear();
  patElecDr03HcalDepth2TowerSumEt_.clear();
  patElecDr03HcalTowerSumEt_.clear();
  patElecDr03TkSumPt_.clear();

  isPassVeto_.clear();
  isPassLoose_.clear();
  isPassMedium_.clear();
  isPassTight_.clear();
  isPassHEEP_.clear();
  isPassVetoNoIso_.clear();
  isPassLooseNoIso_.clear();
  isPassMediumNoIso_.clear();
  isPassTightNoIso_.clear();
  isPassHEEPNoIso_.clear();
  isPassMVAMedium_.clear();
  isPassMVATight_.clear();

  mvaValue_.clear();
  mvaCategory_.clear();

}
