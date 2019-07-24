/* 
-- added CA15 double b-tagger 
-- added ECFs
-- including many files from the external packages of CMSSW (fastjet) to recluster the jet. 
 */
#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/BTauReco/interface/SecondaryVertexTagInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "ExoPieElement/TreeMaker/interface/jetTree.h"
#include <CLHEP/Vector/LorentzVector.h>
#include "TMath.h"
#include "Math/VectorUtil.h"
#include <algorithm>

// for CA15 double b-tagger variables. 
#include "DataFormats/BTauReco/interface/TaggingVariable.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/BTauReco/interface/BoostedDoubleSVTagInfo.h"


// Fastjet 
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/MeasureDefinition.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"

const double DUMMY=-99999.;


typedef math::XYZTLorentzVector LorentzVector;
const math::XYZPoint & position(const reco::Vertex & sv) {return sv.position();}
const math::XYZPoint & position(const reco::VertexCompositePtrCandidate & sv) {return sv.vertex();}

jetTree::jetTree(std::string desc, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(desc, tree),
  isTHINJet_(false),
  isFATJet_(false),
  isADDJet_(false),
  isAK4PuppiJet_(false),
  isAK8PuppiJet_(false),
  isCA15PuppiJet_(false),
  useJECText_(iConfig.getParameter<bool>("useJECText")),
  svTagInfosCstr_(iConfig.getParameter<std::string>("svTagInfosPY")),
  jecUncPayLoadName_(iConfig.getParameter<std::string>(Form("%sjecUncPayLoad",desc.data()))),
  jecNames_(iConfig.getParameter<std::vector<std::string> >(Form("%sjecNames",desc.data()) )), 
  jecUncName_(iConfig.getParameter<std::string>(Form("%sjecUncName",desc.data())) ),	
  jet2012ID_()
{
  
  if (desc.find("THIN")!=std::string::npos)
    isTHINJet_=true;
  if (desc.find("FAT")!=std::string::npos)
    isFATJet_=true;
  if (desc.find("ADD")!=std::string::npos)
    isADDJet_=true; 
  if (desc.find("AK4Puppi")!=std::string::npos)
    isAK4PuppiJet_=true; 
  if (desc.find("AK8Puppi")!=std::string::npos)
    isAK8PuppiJet_=true; 
  if (desc.find("CA15Puppi")!=std::string::npos)
    isCA15PuppiJet_=true; 

  std::cout << desc << std::endl;
  
  
  genjetP4_    = new TClonesArray("TLorentzVector");
  jetP4_       = new TClonesArray("TLorentzVector");
  unCorrJetP4_ = new TClonesArray("TLorentzVector");
  jetPuppiP4_  = new TClonesArray("TLorentzVector");
  jetPuppiSDRawP4_  = new TClonesArray("TLorentzVector");

  SetBranches();


  if(isFATJet_)
    {
      prunedMassJecNames_          = iConfig.getParameter<std::vector<std::string> >(Form("%sprunedMassJecNames",desc.data()));

      if(useJECText_){

	std::vector<JetCorrectorParameters> vPar;

	// pruned mass

	for ( std::vector<std::string>::const_iterator payloadBegin = 
		prunedMassJecNames_.begin(),
		payloadEnd = prunedMassJecNames_.end(), ipayload = payloadBegin; 
	      ipayload != payloadEnd; ++ipayload ) 
	  {
	    JetCorrectorParameters pars(*ipayload);
	    vPar.push_back(pars);
	  }
	prunedjecText_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );


      }

    } // if it's FATjet

  if(isCA15PuppiJet_)
    {
      std::string cmssw_base = getenv("CMSSW_BASE");
      std::string fweight = cmssw_base+"/src/ExoPieElement/TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml";
      mJetBoostedBtaggingMVACalc.initialize("BDT",fweight);

      /* ECF: Starts here */
      jetDefCA = new fastjet::JetDefinition(fastjet::cambridge_algorithm, radius);
      double sdZcut, sdBeta;
      if (radius<1) {
	sdZcut=0.1; sdBeta=0.;
      } else {
	sdZcut=0.15; sdBeta=1.;
      }
      softdrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

      int activeAreaRepeats = 1;
      double ghostArea = 0.01;
      double ghostEtaMax = 7.0;
      activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
      areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);
      ecfnmanager = new ECFNManager();

      /* ECF: Starts here */
    } //  if(isCA15PuppiJet_)


  if(useJECText_)
    {
      
      std::vector<JetCorrectorParameters> vPar;

      for ( std::vector<std::string>::const_iterator payloadBegin = 
	      jecNames_.begin(),
	      payloadEnd = jecNames_.end(), ipayload = payloadBegin; 
	    ipayload != payloadEnd; ++ipayload ) 
	{
	  JetCorrectorParameters pars(*ipayload);
	  vPar.push_back(pars);
	}
      jecText_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );

      jecUncText_ = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecUncName_) );
    }

  


}


jetTree::~jetTree(){

  delete genjetP4_;
  delete jetP4_;
  delete unCorrJetP4_;
  delete jetPuppiP4_;
  delete jetPuppiSDRawP4_;
  
  /* EFC: starts here */
  delete areaDef;
  delete activeArea;
  delete jetDefCA;
  delete softdrop;
  delete tau;
  delete ecfnmanager;
  delete softdrop;
  /* EFC: ends here */

  //delete htt;
  //delete mMCJetCorrector;
  //delete mDataJetCorrector;

}


void
jetTree::Fill(const edm::Event& iEvent, edm::EventSetup const& iSetup){
  Clear();
 
  // // Get the rho collection
  edm::Handle< double > h_rho;
  if(not iEvent.getByToken(rhoForJetToken, h_rho ))
    {
      std::cout<<"FATAL EXCEPTION: in beginging "<<"Following Not Found: "
	       <<"fixedGridRhoFastjetAll" <<std::endl; 
      exit(0);
    }

  jetRho_ = *(h_rho.product());

  // // Get the primary vertex collection                                         
  edm::Handle<reco::VertexCollection>  h_pv;
  if(not iEvent.getByToken(vertexToken,h_pv))
    {
      std::cout<<"FATAL EXCEPTION: in beginging "<<"Following Not Found: "
	       <<"vertexToken"<<std::endl; 
      exit(0);
    }

  if (h_pv->empty()) return; // skip the event if no PV found

  jetNPV_=  h_pv->size();


 
  edm::Handle<pat::JetCollection> JetHandle;
  if(not iEvent.getByToken(jetToken,JetHandle))
    {
      std::cout<<"FATAL EXCEPTION: in beginging "<<"Following Not Found: "
	       <<"jetToken"<<std::endl; 
      exit(0);
    }


  // for getting the L2+L3 correction factor of pruned jet mass
  edm::Handle<pat::JetCollection> JetHandleForPrunedMass;
  pat::JetCollection jetsForPrunedMass;

  if(isFATJet_ && not iEvent.getByToken(prunedMToken,JetHandleForPrunedMass))
    {
      std::cout<<"FATAL EXCEPTION: in beginging "<<"Following Not Found: "
    	       <<"PrunedMassJet"<<std::endl; 
      exit(0);
    }
  else if(isFATJet_ && iEvent.getByToken(prunedMToken,JetHandleForPrunedMass))
    jetsForPrunedMass       = *(JetHandleForPrunedMass.product());



  // for jet energy uncertainty, using global tag
  JetCorrectionUncertainty *jecUnc_=0;
  // fat jet uncertainty does not exist yet
  if(!useJECText_){
    edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
    iSetup.get<JetCorrectionsRecord>().get(jecUncPayLoadName_.data(),JetCorParColl); 
    JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
    jecUnc_ = new JetCorrectionUncertainty(JetCorPar);
  }

  // now start looping over jets
  pat::JetCollection jets(*(JetHandle.product()));
  std::sort(jets.begin(),jets.end(),PtGreater());
  std::vector<pat::Jet>::const_iterator jet =jets.begin();   


  for(;jet!=jets.end();jet++){

    if(jet->pt() < 10.) continue;
  
    nJet_++;
    //Stuff common for all jets.
    

    if (jet->genJet()){

      double DR =0 ;
      double a = ( jet->p4().eta() - jet->genJet()->p4().eta()) * ( jet->p4().eta() - jet->genJet()->p4().eta());
      double b = ( jet->p4().phi() - jet->genJet()->p4().phi()) * ( jet->p4().phi() - jet->genJet()->p4().phi());
      DR = sqrt(a+b);
    
      new( (*genjetP4_)[nJet_-1]) TLorentzVector(
   						 jet->genJet()->p4().px(),
   						 jet->genJet()->p4().py(),
   						 jet->genJet()->p4().pz(),
   						 jet->genJet()->p4().energy()
   						 );
    
      genjetEM_.push_back(jet->genJet()->emEnergy());
      genjetHAD_.push_back(jet->genJet()->hadEnergy());
      genjetINV_.push_back(jet->genJet()->invisibleEnergy());
      genjetAUX_.push_back(jet->genJet()->auxiliaryEnergy());
      matchedDR_.push_back(DR);
      if(false) std::cout<<" jet E = "<<jet->energy()
   			 <<" genjet E = "<<jet->genJet()->energy()
   			 <<" genjet em = "<<jet->genJet()->emEnergy()
   			 <<" genjet had= "<<jet->genJet()->hadEnergy()
   			 <<" genjet inv= "<<jet->genJet()->invisibleEnergy()
   			 <<" genjet aux= "<<jet->genJet()->auxiliaryEnergy()
   			 <<" recojet pho = "<<jet->photonEnergyFraction()*jet->energy()
   			 <<" recojet ele = "<<jet->chargedEmEnergyFraction()*jet->energy()
   			 <<" dr = "<<DR
   			 <<std::endl;
    }
    else{
      genjetEM_.push_back(DUMMY);
      genjetHAD_.push_back(DUMMY);
      genjetINV_.push_back(DUMMY);
      genjetAUX_.push_back(DUMMY);
      matchedDR_.push_back(DUMMY);
      new( (*genjetP4_)[nJet_-1]) TLorentzVector(DUMMY,DUMMY,DUMMY,DUMMY);
    }
  
    jetRawFactor_.push_back(jet->jecFactor("Uncorrected"));
    reco::Candidate::LorentzVector uncorrJet;
    uncorrJet = jet->correctedP4(0);
    new( (*unCorrJetP4_)[nJet_-1]) TLorentzVector(	
						  uncorrJet.px(),
						  uncorrJet.py(),
						  uncorrJet.pz(),
						  uncorrJet.energy()
   							);
   
    jetArea_.push_back(jet->jetArea());
  
    // if reading text files, set jet 4-momentum
    // make correction using jecText files
    if(useJECText_){
      jecText_->setJetEta( uncorrJet.eta() );
      jecText_->setJetPt ( uncorrJet.pt() );
      jecText_->setJetE  ( uncorrJet.energy() );
      jecText_->setJetA  ( jet->jetArea() );
      jecText_->setRho   ( *(h_rho.product()) );
      jecText_->setNPV   ( h_pv->size() );
      Float_t corr_jet = jecText_->getCorrection();

      new( (*jetP4_)[nJet_-1]) TLorentzVector(uncorrJet.px()*corr_jet,
					      uncorrJet.py()*corr_jet,
					      uncorrJet.pz()*corr_jet,
					      uncorrJet.energy()*corr_jet);
      jecUncText_->setJetEta( uncorrJet.eta() );
      jecUncText_->setJetPt( corr_jet * uncorrJet.pt() );
      jetCorrUncUp_.push_back(jecUncText_->getUncertainty(true));

      jecUncText_->setJetEta( uncorrJet.eta() );
      jecUncText_->setJetPt( corr_jet * uncorrJet.pt() );
      jetCorrUncDown_.push_back(jecUncText_->getUncertainty(false));

    }
    else
      new( (*jetP4_)[nJet_-1]) TLorentzVector(jet->p4().px(),
					      jet->p4().py(),
					      jet->p4().pz(),
					      jet->p4().energy());


    // get jet energy scale uncertainty and related input variables
    // fat jet uncertainty does not exist yet, if using database
    if(!useJECText_){

      jecUnc_->setJetEta(jet->eta());
      jecUnc_->setJetPt(jet->pt()); 
      jetCorrUncUp_.push_back(jecUnc_->getUncertainty(true));

      jecUnc_->setJetEta(jet->eta());
      jecUnc_->setJetPt(jet->pt()); 
      jetCorrUncDown_.push_back(jecUnc_->getUncertainty(false));
    }


    jetCharge_.push_back(jet->charge());
    jetPartonFlavor_.push_back(jet->partonFlavour());
    jetHadronFlavor_.push_back(jet->hadronFlavour());


  
    std::map<std::string, bool> Pass = jet2012ID_.LooseJetCut(*jet);
    bool passOrNot = PassAll(Pass); 
    jetPassIDLoose_.push_back(passOrNot);


    std::map<std::string, bool> PassT = jet2012ID_.TightJetCut(*jet);
    bool passOrNotT = PassAll(PassT); 
    jetPassIDTight_.push_back(passOrNotT);


    if(isTHINJet_){
      float jpumva=0.;
      jpumva= jet->userFloat("pileupJetId:fullDiscriminant");
      //std::cout<<" jpumva = "<<jpumva<<std::endl;
      PUJetID_.push_back(jpumva);
          
      // float jpt = jet->pt();
      // float jeta = jet->eta();

      isPUJetIDLoose_.push_back(
				bool(jet->userInt("pileupJetId:fullId") & (1 << 2)));
      isPUJetIDMedium_.push_back(
				 bool(jet->userInt("pileupJetId:fullId") & (1 << 1)));
      isPUJetIDTight_.push_back(
				bool(jet->userInt("pileupJetId:fullId") & (1 << 0)));
    }
        
    jetCEmEF_.push_back(jet->chargedEmEnergyFraction());
    jetCHadEF_.push_back(jet->chargedHadronEnergyFraction());
    jetPhoEF_.push_back(jet->photonEnergyFraction());
    jetNEmEF_.push_back(jet->neutralEmEnergyFraction());
    jetNHadEF_.push_back(jet->neutralHadronEnergyFraction());

    jetEleEF_.push_back(jet->electronEnergyFraction());
    jetMuoEF_.push_back(jet->muonEnergyFraction());
    jetChMuEF_.push_back(jet->chargedMuEnergyFraction());
      
    jetHFHadEF_.push_back(jet->HFHadronEnergyFraction());
    jetHFEMEF_.push_back(jet->HFEMEnergyFraction());
    jetHOEnergy_.push_back(jet->hoEnergy());
    jetHOEF_.push_back(jet->hoEnergyFraction());
      
    if(false) std::cout<<"jetHFHadEF_ = "<<(jet->HFHadronEnergyFraction())
		       <<"  jetHFEMEF_ = "<<(jet->HFEMEnergyFraction())
		       <<"  jetCHHadMultiplicity_ = "<<(jet->chargedHadronMultiplicity())
		       <<"  jetNHadMulplicity_ = "<<(jet->neutralHadronMultiplicity())
		       <<"  jetPhMultiplicity_ = "<<(jet->photonMultiplicity())
		       <<"  jetEleMultiplicity_ = "<<(jet->electronMultiplicity())
		       <<"  jetHFHadMultiplicity_ = "<<(jet->HFHadronMultiplicity())
		       <<"  jetHFEMMultiplicity_ = "<<(jet->HFEMMultiplicity())
		       <<"  jetChMuEF_ = "<<(jet->chargedMuEnergyFraction())
		       <<"  jetNMultiplicity_ = "<<(jet->neutralMultiplicity())
		       <<"  jetHOEnergy_ = "<<(jet->hoEnergy())
		       <<"  jetHOEF_ = "<<(jet->hoEnergyFraction())
		       <<std::endl;

    jetCMulti_.push_back(jet->chargedMultiplicity());      
    jetEleMultiplicity_.push_back(jet->electronMultiplicity());
    jetMuoMultiplicity_.push_back(jet->muonMultiplicity());

    jetCHHadMultiplicity_.push_back(jet->chargedHadronMultiplicity());
    jetPhMultiplicity_.push_back(jet->photonMultiplicity());
    jetNMultiplicity_.push_back(jet->neutralMultiplicity());
    jetNHadMulplicity_.push_back(jet->neutralHadronMultiplicity());
    jetHFHadMultiplicity_.push_back(jet->HFHadronMultiplicity());
    jetHFEMMultiplicity_.push_back(jet->HFEMMultiplicity());

    


    // b-tagging

    jetSSV_.push_back(jet->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
    jetCSV_.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    jetDeepCSV_.push_back(jet->bDiscriminator("deepFlavourJetTags:probb"));
    jetSSVHE_.push_back(jet->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));      
    jetCISVV2_.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetTCHP_.push_back(jet->bDiscriminator("pfTrackCountingHighPurBJetTags"));
    jetTCHE_.push_back(jet->bDiscriminator("pfTrackCountingHighEffBJetTags"));
    jetJP_.push_back(jet->bDiscriminator("pfJetProbabilityBJetTags"));
    jetJBP_.push_back(jet->bDiscriminator("pfJetBProbabilityBJetTags"));


    if(isAK8PuppiJet_){
      jetTau1_.push_back(jet->userFloat("NjettinessAK8Puppi:tau1"));
      jetTau2_.push_back(jet->userFloat("NjettinessAK8Puppi:tau2"));
      jetTau3_.push_back(jet->userFloat("NjettinessAK8Puppi:tau3"));
      jetTau21_.push_back(jet->userFloat("NjettinessAK8Puppi:tau2")/jet->userFloat("NjettinessAK8Puppi:tau1"));
      jetSDmass_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropMass"));
    }


    if(isCA15PuppiJet_){
      
      // For ECFs 
      // reset the ECFs  
      // std::cout<<" now starting ECF"<<std::endl;
      PFatJet *p_jet = new PFatJet(); // choose a better place for this..
      std::vector<float> betas = {0.5,1.,2.,4.};
      std::vector<int> Ns = {1,2,3,4};
      std::vector<int> orders = {1,2,3};
      for (unsigned int iB=0; iB!=4; ++iB) {
        for (auto N : Ns) {
          for (auto o : orders) {
	    p_jet->set_ecf(o,N,iB,-1);
          }
        }
      }
      std::vector<edm::Ptr<reco::Candidate>> constituentPtrs = jet->getJetConstituents();
      /* commented by raman
      if (pfcands!=0) { // associate to pf cands in tree
	const std::map<const reco::Candidate*, UShort_t> &pfmap = pfcands->get_map();
	jet->constituents = new std::vector<UShort_t>();
	std::vector<UShort_t> *constituents = jet->constituents;
	
	for (auto ptr : constituentPtrs) {
	  const reco::Candidate *constituent = ptr.get();
	  
	  auto result_ = pfmap.find(constituent); // check if we know about this pf cand
	  if (result_ == pfmap.end()) {
	    PError("PandaProdNtupler::FatJetFiller",TString::Format("could not PF [%s] ...\n",treename.Data()));
	  } else {
	    constituents->push_back(result_->second);
	  } 
	} // loop through constituents from input
      } // pfcands!=0 
      
      */
      // if (!minimal && data->size()<2)   commented by raman
      {
	// calculate ECFs, groomed tauN
	typedef std::vector<fastjet::PseudoJet> VPseudoJet;
	VPseudoJet vjet;
	for (auto ptr : constituentPtrs) { 
	  // create vector of PseudoJets
	  const reco::Candidate *constituent = ptr.get();
	  if (constituent->pt()<0.01) 
	    continue;
	  vjet.emplace_back(constituent->px(),constituent->py(),constituent->pz(),constituent->energy());
	}
	
	fastjet::ClusterSequenceArea seq(vjet, *jetDefCA, *areaDef); 
	VPseudoJet alljets = fastjet::sorted_by_pt(seq.inclusive_jets(0.1));
	
	if (alljets.size()>0){
	  fastjet::PseudoJet *leadingJet = &(alljets[0]);
	  fastjet::PseudoJet sdJet = (*softdrop)(*leadingJet);
	  // get and filter constituents of groomed jet
	  VPseudoJet sdconsts = fastjet::sorted_by_pt(sdJet.constituents());
	  int nFilter = TMath::Min(100,(int)sdconsts.size());
	  VPseudoJet sdconstsFiltered(sdconsts.begin(),sdconsts.begin()+nFilter);
	  // calculate ECFs
	  for (unsigned int iB=0; iB!=4; ++iB) {
	    calcECFN(betas[iB],sdconstsFiltered,ecfnmanager); // calculate for all Ns and os
	    for (auto N : Ns) {
	      for (auto o : orders) {
		float x = ecfnmanager->ecfns[TString::Format("%i_%i",N,o)];
		int r = p_jet->set_ecf(o,N,iB,x);
		/*std::cout<<"r,o,N,beta,x =  "
			 <<" "<<r
			 <<" "<<o
			 <<" "<<N
			 <<" "<<betas[iB]
			 <<" "<<x
			 <<std::endl;*/
		//std::cout<< " using the get function = "<< p_jet->get_ecf(o,N,iB)<<std::endl;

		if (r) {
		  std::cout<<"FatJetFiller eoor "<<std::endl;
		}
	      } // o loop
	    } // N loop
	  } // beta loop
	  
	}/* // by raman
	    jet->tau3SD = tau->getTau(3,sdconsts);
	  jet->tau2SD = tau->getTau(2,sdconsts);
	  jet->tau1SD = tau->getTau(1,sdconsts);
	  
	  // HTT
	  fastjet::PseudoJet httJet = htt->result(*leadingJet);
	  if (httJet!=0) {
	    fastjet::HEPTopTaggerV2Structure *s = 
	      (fastjet::HEPTopTaggerV2Structure*)httJet.structure_non_const_ptr();
	    jet->htt_mass = s->top_mass();
	    jet->htt_frec = s->fRec();
	  }
	  
	} else {
	  PError("PandaProd::Ntupler::FatJetFiller","Jet could not be clustered");
	}
*/	
      } // if not minimal and fewer than 2 
      // End of ECFs computation 
      
      // std::cout<< " using the get function outside  = "<< p_jet->get_ecf(2,3,1)<<"   "<<p_jet->get_ecf(1,2,1)<<std::endl;
      const reco::TaggingVariableList vars = jet->tagInfoBoostedDoubleSV()->taggingVariables();
      float z_ratio_                       = vars.get(reco::btau::z_ratio);
      float trackSipdSig_3_                = vars.get(reco::btau::trackSip3dSig_3);
      float trackSipdSig_2_                = vars.get(reco::btau::trackSip3dSig_2);
      float trackSipdSig_1_                = vars.get(reco::btau::trackSip3dSig_1);
      float trackSipdSig_0_                = vars.get(reco::btau::trackSip3dSig_0);
      float trackSipdSig_1_0_              = vars.get(reco::btau::tau2_trackSip3dSig_0);
      float trackSipdSig_0_0_              = vars.get(reco::btau::tau1_trackSip3dSig_0);
      float trackSipdSig_1_1_              = vars.get(reco::btau::tau2_trackSip3dSig_1);
      float trackSipdSig_0_1_              = vars.get(reco::btau::tau1_trackSip3dSig_1);
      float trackSip2dSigAboveCharm_0_     = vars.get(reco::btau::trackSip2dSigAboveCharm);
      float trackSip2dSigAboveBottom_0_    = vars.get(reco::btau::trackSip2dSigAboveBottom_0);
      float trackSip2dSigAboveBottom_1_    = vars.get(reco::btau::trackSip2dSigAboveBottom_1);
      float tau1_trackEtaRel_0_            = vars.get(reco::btau::tau2_trackEtaRel_0);
      float tau1_trackEtaRel_1_            = vars.get(reco::btau::tau2_trackEtaRel_1);
      float tau1_trackEtaRel_2_            = vars.get(reco::btau::tau2_trackEtaRel_2);
      float tau0_trackEtaRel_0_            = vars.get(reco::btau::tau1_trackEtaRel_0);
      float tau0_trackEtaRel_1_            = vars.get(reco::btau::tau1_trackEtaRel_1);
      float tau0_trackEtaRel_2_            = vars.get(reco::btau::tau1_trackEtaRel_2);
      float tau_vertexMass_0_              = vars.get(reco::btau::tau1_vertexMass);
      float tau_vertexEnergyRatio_0_       = vars.get(reco::btau::tau1_vertexEnergyRatio);
      float tau_vertexDeltaR_0_            = vars.get(reco::btau::tau1_vertexDeltaR);
      float tau_flightDistance2dSig_0_     = vars.get(reco::btau::tau1_flightDistance2dSig);
      float tau_vertexMass_1_              = vars.get(reco::btau::tau2_vertexMass);
      float tau_vertexEnergyRatio_1_       = vars.get(reco::btau::tau2_vertexEnergyRatio);
      float tau_flightDistance2dSig_1_     = vars.get(reco::btau::tau2_flightDistance2dSig);
      float jetNTracks_                    = vars.get(reco::btau::jetNTracks);
      float nSV_                           = vars.get(reco::btau::jetNSecondaryVertices);
      float massPruned_                    = jet->p4().mass(); //jet->m;
      float flavour_                       = -1;   //j.partonFlavor();   // they're spectator variables
      float nbHadrons_                     = -1; //j.hadronFlavor(); // 
      float ptPruned_                      = jet->p4().pt();
      float etaPruned_                     = jet->p4().eta();
      
      std::vector<float> subjetSDCSV_puppi_tmp;
      subjetSDCSV_puppi_tmp.clear();
      auto const & sdSubjetsPuppi = jet->subjets("SoftDrop");
      float SubJet_csv_ = -1.0;
      for ( auto const & it : sdSubjetsPuppi ) {
	subjetSDCSV_puppi_tmp.push_back(it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	SubJet_csv_                    = *(std::max_element(subjetSDCSV_puppi_tmp.begin(), subjetSDCSV_puppi_tmp.end()));
	if ((SubJet_csv_ < -1) || (SubJet_csv_ > 1)) SubJet_csv_ = -1;
      }
      double double_CA15 = mJetBoostedBtaggingMVACalc.mvaValue(
							       massPruned_,
							       flavour_,
							       nbHadrons_,
							       ptPruned_,
							       etaPruned_,
							       SubJet_csv_,
							       z_ratio_,
							       trackSipdSig_3_,
							       trackSipdSig_2_,
							       trackSipdSig_1_,
							       trackSipdSig_0_,
							       trackSipdSig_1_0_,
							       trackSipdSig_0_0_,
							       trackSipdSig_1_1_,
							       trackSipdSig_0_1_,
							       trackSip2dSigAboveCharm_0_,
							       trackSip2dSigAboveBottom_0_,
							       trackSip2dSigAboveBottom_1_,
							       tau0_trackEtaRel_0_,
							       tau0_trackEtaRel_1_,
							       tau0_trackEtaRel_2_,
							       tau1_trackEtaRel_0_,
							       tau1_trackEtaRel_1_,
							       tau1_trackEtaRel_2_,
							       tau_vertexMass_0_,
							       tau_vertexEnergyRatio_0_,
							       tau_vertexDeltaR_0_,
							       tau_flightDistance2dSig_0_,
							       tau_vertexMass_1_,
							       tau_vertexEnergyRatio_1_,
							       tau_flightDistance2dSig_1_,
							       jetNTracks_,
							       nSV_,
							       false
							       );
      
      // std::cout<<" double_CA15 = "<<double_CA15 << std::endl;
      ca15_doublebtag.push_back(double_CA15);
      ECF_2_3_10.push_back(p_jet->get_ecf(2,3,1)) ;
      ECF_1_2_10.push_back(p_jet->get_ecf(1,2,1));
      
      jetTau1_.push_back(jet->userFloat("NjettinessCA15Puppi:tau1"));
      jetTau2_.push_back(jet->userFloat("NjettinessCA15Puppi:tau2"));
      jetTau3_.push_back(jet->userFloat("NjettinessCA15Puppi:tau3"));
      jetTau21_.push_back(jet->userFloat("NjettinessCA15Puppi:tau2")/jet->userFloat("NjettinessCA15Puppi:tau1"));
      jetSDmass_.push_back(jet->userFloat("ca15PFJetsPuppiSoftDropMass"));
      
      
    }//     if(isCA15PuppiJet_){
    
    
    if(isFATJet_){


      jetTau1_.push_back(jet->userFloat("NjettinessAK8:tau1"));
      jetTau2_.push_back(jet->userFloat("NjettinessAK8:tau2"));
      jetTau3_.push_back(jet->userFloat("NjettinessAK8:tau3"));
      jetTau21_.push_back(jet->userFloat("NjettinessAK8:tau2")/jet->userFloat("NjettinessAK8:tau1"));
      

      //Puppi related information
      jetPuppiTau1_.push_back(jet->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));
      jetPuppiTau2_.push_back(jet->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
      jetPuppiTau3_.push_back(jet->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3"));


      TLorentzVector temp_puppi;
      temp_puppi.SetPtEtaPhiM(jet->userFloat("ak8PFJetsPuppiValueMap:pt"),
			      jet->userFloat("ak8PFJetsPuppiValueMap:eta"),
			      jet->userFloat("ak8PFJetsPuppiValueMap:phi"),
			      jet->userFloat("ak8PFJetsPuppiValueMap:mass"));
			      

      new( (*jetPuppiP4_)[nJet_-1]) TLorentzVector(temp_puppi);



      unsigned int nSubSoftDropjets_puppi=0;	
  	
      std::vector<int>   subjetSDFatJetIndex_puppi;
      std::vector<float> subjetSDPx_puppi; 
      std::vector<float> subjetSDPy_puppi; 
      std::vector<float> subjetSDPz_puppi; 
      std::vector<float> subjetSDE_puppi; 	
      std::vector<float> subjetSDCSV_puppi; 	

      subjetSDFatJetIndex_puppi.clear();
      subjetSDPx_puppi.clear();
      subjetSDPy_puppi.clear();
      subjetSDPz_puppi.clear();
      subjetSDE_puppi.clear();
      subjetSDCSV_puppi.clear();

      TLorentzVector puppi_softdrop(0,0,0,0);
      TLorentzVector puppi_softdrop_raw(0,0,0,0);
      auto const & sdSubjetsPuppi = jet->subjets("SoftDropPuppi");
      for ( auto const & it : sdSubjetsPuppi ) {
	nSubSoftDropjets_puppi++;

	subjetSDFatJetIndex_puppi.push_back(nJet_-1);
	subjetSDPx_puppi.push_back(it->px());
	subjetSDPy_puppi.push_back(it->py());
	subjetSDPz_puppi.push_back(it->pz());
	subjetSDE_puppi.push_back(it->energy());	
	subjetSDCSV_puppi.push_back(it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));	

	puppi_softdrop += TLorentzVector(it->px(),
					 it->py(),
					 it->pz(),
					 it->energy());

	puppi_softdrop_raw += TLorentzVector(it->correctedP4(0).px(),
					     it->correctedP4(0).py(),
					     it->correctedP4(0).pz(),
					     it->correctedP4(0).energy());
      } //subjet loop

      new( (*jetPuppiSDRawP4_)[nJet_-1]) TLorentzVector(puppi_softdrop_raw);


      //      using a different way to get corrected pruned/softdrop mass
      // if reading global tag
      float corr=-1;

      if(!useJECText_){	
	// pruned mass: CHS
	std::vector<pat::Jet>::const_iterator jetForPrunedMass = 
	  find_if(jetsForPrunedMass.begin(),
		  jetsForPrunedMass.end(),
		  [&jet](const pat::Jet& item)->bool{return fabs(jet->correctedP4(0).pt()-item.correctedP4(0).pt())<1e-3;});	
    
	if(jetForPrunedMass!=jetsForPrunedMass.end())
	  corr = jetForPrunedMass->pt()/jetForPrunedMass->correctedP4(0).pt();
      }

      else if(useJECText_){
	// pruned mass: CHS
	prunedjecText_->setJetEta( uncorrJet.eta() );
	prunedjecText_->setJetPt ( uncorrJet.pt() );
	prunedjecText_->setJetE  ( uncorrJet.energy() );
	prunedjecText_->setJetA  ( jet->jetArea() );
	prunedjecText_->setRho   ( *(h_rho.product()) );
	prunedjecText_->setNPV   ( h_pv->size() );
	corr = prunedjecText_->getCorrection();

      }

      if(corr<0)
	jetPRmassL2L3Corr_.push_back(DUMMY);
      else
	jetPRmassL2L3Corr_.push_back(corr*jet->userFloat("ak8PFJetsCHSPrunedMass"));


      jetSDmass_.push_back(jet->userFloat("ak8PFJetsCHSSoftDropMass"));
      jetPRmass_.push_back(jet->userFloat("ak8PFJetsCHSPrunedMass"));
 


      if(nSubSoftDropjets_puppi==0)
	{
	  subjetSDFatJetIndex_puppi.push_back(DUMMY);
	  subjetSDPx_puppi.push_back(DUMMY);
	  subjetSDPy_puppi.push_back(DUMMY);
	  subjetSDPz_puppi.push_back(DUMMY);
	  subjetSDE_puppi.push_back(DUMMY);	
	  subjetSDCSV_puppi.push_back(DUMMY);	
	  jetPuppiSDmass_.push_back(DUMMY);
	}
      else
	{
	  jetPuppiSDmass_.push_back(puppi_softdrop_raw.M());
	}

           
      nSubSDPuppiJet_.push_back(nSubSoftDropjets_puppi); 
      subjetSDPuppiFatJetIndex_.push_back(subjetSDFatJetIndex_puppi);
      subjetSDPuppiPx_.push_back(subjetSDPx_puppi);
      subjetSDPuppiPy_.push_back(subjetSDPy_puppi);
      subjetSDPuppiPz_.push_back(subjetSDPz_puppi);
      subjetSDPuppiE_.push_back(subjetSDE_puppi);
      subjetSDPuppiCSV_.push_back(subjetSDCSV_puppi);



    } // only for AK8CHS jets

  
    // if this is a AK8 jet
    if(!isTHINJet_ && !isAK4PuppiJet_){
      //HBB tagger
      if(isCA15PuppiJet_)
	jet_DoubleSV_.push_back(jet->bDiscriminator("pfBoostedDoubleSecondaryVertexCA15BJetTags"));
      else
	jet_DoubleSV_.push_back(jet->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
        
      std::vector<float> jet_SVMass_float;
  	  
      jet_SVMass_float.clear();
      unsigned int nSV=0;  

      if(jet->hasTagInfo(svTagInfosCstr_.data()))
	{
	  const reco::CandSecondaryVertexTagInfo *candSVTagInfo = jet->tagInfoCandSecondaryVertex("pfInclusiveSecondaryVertexFinder");
	  nSV = candSVTagInfo->nVertices();                       
  	  
	  for(unsigned int n_2ndvtx=0;n_2ndvtx< nSV;n_2ndvtx++)
	    jet_SVMass_float.push_back(candSVTagInfo->secondaryVertex(n_2ndvtx).p4().mass());
  	  
	} // if there is tagging information
      if(nSV==0)
	jet_SVMass_float.push_back(DUMMY);	    
      jet_nSV_.push_back(nSV);  
      jet_SVMass_.push_back(jet_SVMass_float);
    
    } /// if its ADDJET or FATjet or AK8PuppiJet or CA15PuppiJet
 
    // softdrop subjets
    if(isFATJet_ || isAK8PuppiJet_ || isCA15PuppiJet_){
      pat::Jet const *jetptr = &*jet;  
      
      //	std::cout<<" working before SoftDrop "<<std::endl;
  	
      auto wSubjets = jetptr->subjets("SoftDrop");
      //	std::cout<<" working after SoftDrop "<<std::endl;
      int nSubSoftDropjets=0;	
  	
      std::vector<int>   subjetSDFatJetIndex;
      std::vector<float> subjetSDPx; 
      std::vector<float> subjetSDPy; 
      std::vector<float> subjetSDPz; 
      std::vector<float> subjetSDE; 	
      std::vector<float> subjetSDRawFactor; 
      std::vector<int>   subjetSDCharge;
      std::vector<int>   subjetSDPartonFlavor;
      std::vector<int>   subjetSDHadronFlavor;
      std::vector<float> subjetSDCSV; 

      subjetSDFatJetIndex.clear();
      subjetSDPx.clear();
      subjetSDPy.clear();
      subjetSDPz.clear();
      subjetSDE.clear();
      subjetSDRawFactor.clear();
      subjetSDCharge.clear();
      subjetSDPartonFlavor.clear();
      subjetSDHadronFlavor.clear();
      subjetSDCSV.clear(); 

      float genjet_softdropmass=DUMMY;
      TLorentzVector genjet_softdrop_l4;
      genjet_softdrop_l4.SetPxPyPzE(0,0,0,0);

      for ( auto const & iw : wSubjets ) 
	{

	  nSubSoftDropjets++;

	  // build genjet softdrop mass
	  if (iw->genJet()){
	    genjet_softdrop_l4 +=TLorentzVector(
						iw->genJet()->p4().px(),
						iw->genJet()->p4().py(),
						iw->genJet()->p4().pz(),
						iw->genJet()->p4().energy()
						);
	  }
  	      
	  subjetSDFatJetIndex.push_back(nJet_-1);

	  subjetSDPx.push_back(iw->px());
	  subjetSDPy.push_back(iw->py());
	  subjetSDPz.push_back(iw->pz());
	  subjetSDE.push_back(iw->energy());	
	  subjetSDRawFactor.push_back(iw->jecFactor("Uncorrected"));
	  subjetSDCharge.push_back(iw->charge());
	  subjetSDPartonFlavor.push_back(iw->partonFlavour());
	  subjetSDHadronFlavor.push_back(iw->hadronFlavour());	      
	  subjetSDCSV.push_back(iw->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));   

	}//subjet loop
      if(nSubSoftDropjets==0)
	{
	  subjetSDFatJetIndex.push_back(DUMMY);
	  subjetSDPx.push_back(DUMMY);
	  subjetSDPy.push_back(DUMMY);
	  subjetSDPz.push_back(DUMMY);
	  subjetSDE.push_back(DUMMY);	
	  subjetSDRawFactor.push_back(DUMMY);
	  subjetSDCharge.push_back(DUMMY);
	  subjetSDPartonFlavor.push_back(DUMMY);
	  subjetSDHadronFlavor.push_back(DUMMY);	      
	  subjetSDCSV.push_back(DUMMY);   
	}
      
      if(nSubSoftDropjets==0 || genjet_softdrop_l4.E()<1e-6)
	genjet_softdropmass = DUMMY;
      else
	genjet_softdropmass = genjet_softdrop_l4.M();

      jetGenSDmass_.push_back(genjet_softdropmass);
      nSubSDJet_.push_back(nSubSoftDropjets); 
      subjetSDFatJetIndex_.push_back(subjetSDFatJetIndex);
      subjetSDPx_.push_back(subjetSDPx);
      subjetSDPy_.push_back(subjetSDPy);
      subjetSDPz_.push_back(subjetSDPz);
      subjetSDE_.push_back(subjetSDE);
      subjetSDRawFactor_.push_back(subjetSDRawFactor);
      subjetSDCharge_.push_back(subjetSDCharge);
      subjetSDPartonFlavor_.push_back(subjetSDPartonFlavor);
      subjetSDHadronFlavor_.push_back(subjetSDHadronFlavor);	      
      subjetSDCSV_.push_back(subjetSDCSV); 
     
  	  

    }//if is Fat jet  or AK8Puppijet or CA15Puppijet


  }//jet loop
    

  // fat jet uncertainty does not exist yet
  if(!useJECText_)
    delete jecUnc_;


}



void
jetTree::SetBranches(){
  
  AddBranch(&nJet_,   "nJet");
  AddBranch(&jetP4_,       "jetP4");

  if(!isADDJet_){
  AddBranch(&jetRho_, "jetRho");
  AddBranch(&jetNPV_, "jetNPV");

  AddBranch(&genjetP4_,   "genjetP4");
  AddBranch(&genjetEM_ ,  "genjetEM");
  AddBranch(&genjetHAD_ , "genjetHAD");
  AddBranch(&genjetINV_ , "genjetINV");
  AddBranch(&genjetAUX_ , "genjetAUX");
  AddBranch(&matchedDR_ , "matchedDR");

  AddBranch(&jetRawFactor_, "jetRawFactor");
  AddBranch(&unCorrJetP4_, "unCorrJetP4");

  AddBranch(&jetArea_,        "jetArea");
  AddBranch(&jetCorrUncUp_,   "jetCorrUncUp");
  AddBranch(&jetCorrUncDown_, "jetCorrUncDown");

  AddBranch(&jetCharge_,       "jetCharge");
  AddBranch(&jetPartonFlavor_, "jetPartonFlavor");
  AddBranch(&jetHadronFlavor_, "jetHadronFlavor");
  AddBranch(&jetPassIDLoose_,  "jetPassIDLoose");
  AddBranch(&jetPassIDTight_,  "jetPassIDTight");

  AddBranch(&jetCEmEF_,  "jetCEmEF");
  AddBranch(&jetCHadEF_, "jetCHadEF");
  AddBranch(&jetPhoEF_,  "jetPhoEF");
  AddBranch(&jetNEmEF_,  "jetNEmEF");
  AddBranch(&jetNHadEF_, "jetNHadEF");
  AddBranch(&jetEleEF_,  "jetEleEF");
  AddBranch(&jetMuoEF_,  "jetMuoEF");

  AddBranch(&jetCMulti_, "jetCMulti");
  AddBranch(&jetEleMultiplicity_,"jetEleMulti");
  AddBranch(&jetMuoMultiplicity_,"jetMuoMulti");
  
  AddBranch(&jetSSV_,   "jetSSV");
  AddBranch(&jetCSV_,   "jetCSV");        
  AddBranch(&jetSSVHE_, "jetSSVHE");
  AddBranch(&jetCISVV2_,"jetCISVV2");
  AddBranch(&jetTCHP_,  "jetTCHP");
  AddBranch(&jetTCHE_,  "jetTCHE");
  AddBranch(&jetJP_,    "jetJP");
  AddBranch(&jetJBP_,   "jetJBP");


  }

  if(isTHINJet_){
    AddBranch(&PUJetID_,   "PUJetID");
    AddBranch(&jetDeepCSV_,  "jetDeepCSV");
    AddBranch(&isPUJetIDLoose_,  "isPUJetIDLoose");
    AddBranch(&isPUJetIDMedium_, "isPUJetIDMedium");
    AddBranch(&isPUJetIDTight_,  "isPUJetIDTight");
  }

  if(isFATJet_ || isAK8PuppiJet_ || isCA15PuppiJet_){

    AddBranch(&jetTau1_,  "jetTau1");
    AddBranch(&jetTau2_,  "jetTau2");
    AddBranch(&jetTau3_,  "jetTau3");
    AddBranch(&jetTau21_, "jetTau21");
    AddBranch(&jetSDmass_,         "jetSDmass");


    // subjet information
    AddBranch(&jetGenSDmass_,         "jetGenSDmass");
    AddBranch(&nSubSDJet_,            "nSubSDJet");
    AddBranch(&subjetSDFatJetIndex_,  "subjetSDFatJetIndex");
    AddBranch(&subjetSDPx_,           "subjetSDPx");     
    AddBranch(&subjetSDPy_,           "subjetSDPy");     
    AddBranch(&subjetSDPz_,           "subjetSDPz");     
    AddBranch(&subjetSDE_,            "subjetSDE");     
    AddBranch(&subjetSDRawFactor_,    "subjetSDRawFactor");     
    AddBranch(&subjetSDPartonFlavor_, "subjetSDPartonFlavor");
    AddBranch(&subjetSDHadronFlavor_, "subjetSDHadronFlavor");
    AddBranch(&subjetSDCSV_,          "subjetSDCSV");     

    if(isCA15PuppiJet_){
      AddBranch(&ca15_doublebtag, "_doublebtag");
      AddBranch(&ECF_2_3_10, "ECF_2_3_10");
      AddBranch(&ECF_1_2_10, "ECF_1_2_10");
    }
    
  }
  
  if(isFATJet_){

    
    AddBranch(&jetPRmass_,         "jetPRmass");
    AddBranch(&jetPRmassL2L3Corr_, "jetPRmassL2L3Corr");


    // puppi information
    AddBranch(&jetPuppiTau1_,   "jetPuppiTau1");
    AddBranch(&jetPuppiTau2_,   "jetPuppiTau2");
    AddBranch(&jetPuppiTau3_,   "jetPuppiTau3");

    AddBranch(&jetPuppiSDmass_,         "jetPuppiSDmass");

    AddBranch(&jetPuppiP4_, "jetPuppiP4");
    AddBranch(&jetPuppiSDRawP4_, "jetPuppiSDRawP4");

    AddBranch(&nSubSDPuppiJet_,           "nSubSDPuppiJet");
    AddBranch(&subjetSDPuppiFatJetIndex_, "subjetSDPuppiFatJetIndex");
    AddBranch(&subjetSDPuppiPx_,          "subjetSDPuppiPx");     
    AddBranch(&subjetSDPuppiPy_,          "subjetSDPuppiPy");     
    AddBranch(&subjetSDPuppiPz_,          "subjetSDPuppiPz");     
    AddBranch(&subjetSDPuppiE_,           "subjetSDPuppiE");     
    AddBranch(&subjetSDPuppiCSV_,           "subjetSDPuppiCSV");     

  }
  
  if(!isTHINJet_ && !isAK4PuppiJet_)
    {

      AddBranch(&jet_DoubleSV_,"jet_DoubleSV");
      AddBranch(&jet_nSV_,     "jet_nSV");
      AddBranch(&jet_SVMass_,  "jet_SVMass");
    }



}


void
jetTree::Clear(){


  nJet_ =0;
  jetRho_=0;
  jetNPV_=0;

  genjetP4_->Clear();
  genjetEM_.clear();
  genjetHAD_.clear();
  genjetINV_.clear();
  genjetAUX_.clear();
  matchedDR_.clear();

  jetRawFactor_.clear();

  jetP4_->Clear();  
  unCorrJetP4_->Clear();

  jetArea_.clear();
  jetCorrUncUp_.clear();
  jetCorrUncDown_.clear();
  jetCharge_.clear();
  jetPartonFlavor_.clear();
  jetHadronFlavor_.clear();
  jetPassIDLoose_.clear();
  jetPassIDTight_.clear();
  PUJetID_.clear();
  isPUJetIDLoose_.clear();
  isPUJetIDMedium_.clear();
  isPUJetIDTight_.clear();

  //Energy Fraction and Multiplicity 

  jetCEmEF_.clear();
  jetCHadEF_.clear();
  jetPhoEF_.clear();
  jetNEmEF_.clear();
  jetNHadEF_.clear();

  jetMuoEF_.clear();
  jetEleEF_.clear();
  jetChMuEF_.clear();

  jetHFHadEF_.clear();
  jetHFEMEF_.clear();
  jetHOEnergy_.clear();
  jetHOEF_.clear();


  jetCMulti_.clear();
  jetEleMultiplicity_.clear();
  jetMuoMultiplicity_.clear();
  jetCHHadMultiplicity_.clear();
  jetPhMultiplicity_.clear();
  jetNMultiplicity_.clear();
  jetNHadMulplicity_.clear();
  jetHFHadMultiplicity_.clear();
  jetHFEMMultiplicity_.clear();


  // btag information
  jetSSV_.clear();
  jetCSV_.clear();
  jetDeepCSV_.clear();
  jetSSVHE_.clear();
  jetCISVV2_.clear();
  jetTCHP_.clear();
  jetTCHE_.clear();
  jetJP_.clear();
  jetJBP_.clear();


  jetTau1_.clear();
  jetTau2_.clear();
  jetTau3_.clear();
  jetTau21_.clear();


  //ak8jet mass
 
  jetSDmass_.clear(); 
  jetPRmass_.clear();  
  jetPRmassL2L3Corr_.clear();

  // puppi related stuff
  
  jetPuppiTau1_.clear();
  jetPuppiTau2_.clear();
  jetPuppiTau3_.clear();
  jetPuppiSDmass_.clear();

  jetPuppiP4_->Clear();
  jetPuppiSDRawP4_->Clear();
  nSubSDPuppiJet_.clear();
  subjetSDPuppiFatJetIndex_.clear(); 
  subjetSDPuppiPx_.clear();
  subjetSDPuppiPy_.clear();
  subjetSDPuppiPz_.clear();
  subjetSDPuppiE_.clear();
  subjetSDPuppiCSV_.clear();

  // CA15 and ECFs 
  ca15_doublebtag.clear();
  ECF_2_3_10.clear();
  ECF_1_2_10.clear();
  

  //jet  Hbb tagger for fat and add jet

  jet_DoubleSV_.clear();


  //jet secondary vtx

  jet_nSV_.clear();
  jet_SVMass_.clear();
  jet_SVEnergyRatio_.clear();




  // subjet of jets
  jetGenSDmass_.clear();
  nSubSDJet_.clear();
  subjetSDPx_.clear();
  subjetSDPy_.clear();
  subjetSDPz_.clear();
  subjetSDE_.clear();
  subjetSDRawFactor_.clear();
  subjetSDCharge_.clear();
  subjetSDFatJetIndex_.clear();
  subjetSDPartonFlavor_.clear();
  subjetSDHadronFlavor_.clear();
  subjetSDCSV_.clear();        



}
