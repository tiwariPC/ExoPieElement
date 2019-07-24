#include "ExoPieElement/TreeMaker/interface/genInfoTree.h"
#include <bitset>
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "TLorentzVector.h"
const double DUMMY=-99999.;



genInfoTree::genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(name,tree),
  MAXNGENPAR_(iConfig.getParameter<unsigned int>("maxNumGenPar")),
  applyStatusSelection_(iConfig.getParameter<bool>("applyStatusSelection")),
  applyPromptSelection_(iConfig.getParameter<bool>("applyPromptSelection")),
  saveLHEWeights_(iConfig.getParameter<bool>("saveLHEWeights")),
  saveGenJets_(iConfig.getParameter<bool>("saveGenJets")),
  saveGenJetSub_(iConfig.getParameter<bool>("saveGenJetSub"))
{
  genParP4_ =   new TClonesArray("TLorentzVector");
  ak4GenJetP4_ =   new TClonesArray("TLorentzVector");
  ak8GenJetP4_ =   new TClonesArray("TLorentzVector");
  
  SetBranches();

  if(saveGenJets_){

    const double radius=0.8, sdZcut=0.1, sdBeta=0.;

    jetDefAKT = new fastjet::JetDefinition(fastjet::antikt_algorithm, radius);
    softdrop = new fastjet::contrib::SoftDrop(sdBeta,sdZcut,radius);

    fjtau1 = new fastjet::contrib::Nsubjettiness(1, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1,radius));
    fjtau2 = new fastjet::contrib::Nsubjettiness(2, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1,radius));
    fjtau3 = new fastjet::contrib::Nsubjettiness(3, fastjet::contrib::OnePass_KT_Axes(), fastjet::contrib::NormalizedMeasure(1,radius));
    

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    activeArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
    areaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*activeArea);

  } // if(saveGenJets_) ends here

}


genInfoTree::~genInfoTree()
{
  delete genParP4_;
  delete ak4GenJetP4_;
  delete ak8GenJetP4_;

}


//
// member functions
//

void 
genInfoTree::GetRunInfo(const edm::Run& iRun)
{
  edm::Handle<LHERunInfoProduct> run;
  typedef std::vector<LHERunInfoProduct::Header>::const_iterator headers_const_iterator;

  if(iRun.getByToken( lheRunToken, run )){

    LHERunInfoProduct myLHERunInfoProduct = *(run.product());
    for (headers_const_iterator iter=myLHERunInfoProduct.headers_begin(); iter!=myLHERunInfoProduct.headers_end(); iter++){
      // std::cout << iter->tag() << std::endl;                                                                                                                                                              
      std::vector<std::string> lines = iter->lines();
      for (unsigned int iLine = 0; iLine<lines.size(); iLine++) {
	std::string thisline = lines.at(iLine);
	// only print out lines with weight information, otherwise the log file will be too big
        if(thisline.find("weight id")!=std::string::npos ||
           thisline.find("weightgroup")!=std::string::npos)
	  std::cout << thisline;
      }

    }

  }

}



// ------------ method called to for each event  ------------
void
genInfoTree::Fill(const edm::Event& iEvent)
{
  Clear();
  if(iEvent.isRealData())return;

  using namespace edm;
  edm::Handle<reco::GenParticleCollection> genParticleHandle;
  if(not iEvent.getByToken(genParticleToken, genParticleHandle))
    {
      std::cout<<
	"GenAnalyzer: Generator Level Information not found\n"
	       <<std::endl;
    }

  edm::Handle<GenEventInfoProduct>    genEventInfoHandle;

  if (iEvent.getByToken(genEventToken, genEventInfoHandle)) {
    if (genEventInfoHandle->hasBinningValues())
      ptHat_ = genEventInfoHandle->binningValues()[0];
      
    mcWeight_ = genEventInfoHandle->weight();

    if (genEventInfoHandle->pdf()) {
      pdf_.push_back(genEventInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->x.first);     // x value of parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->x.second);    // x value of parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->scalePDF);    // scale of the hard interaction
    }
  }

  // add HT information
  edm::Handle<LHEEventProduct> evt;
  

  if(iEvent.getByToken( lheEventToken, evt )){

    // nominal LHE weight for this event
    originalLHEweight_ = evt->originalXWGTUP(); 

    // get PDF and scale weights
    if(saveLHEWeights_){
      for (unsigned int i=0; i< evt->weights().size(); i++) {      
	float tempMCWeight = genEventInfoHandle->weight() > 0? 1: -1;
	float sysLHEweight = tempMCWeight* evt->weights()[i].wgt/evt->originalXWGTUP();
	pdfscaleSysWeights_.push_back( sysLHEweight );
      }
    }


    HT_=0;
    const lhef::HEPEUP hepeup_ = evt->hepeup();

    const int nup_ = hepeup_.NUP; 
    const std::vector<int> idup_ = hepeup_.IDUP;
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP;
    const std::vector<int> istup_ = hepeup_.ISTUP;

    for ( unsigned int icount = 0 ; icount < (unsigned int)nup_; icount++ ) {

      int PID    = idup_[icount];
      int status = istup_[icount];
      double px = (pup_[icount])[0];
      double py = (pup_[icount])[1];
            
      if(status!=1)continue;
      // if it's not a gluon or quark
      if(!(abs(PID)==21 || (abs(PID)<6 && abs(PID)>0)))continue;
      HT_ += sqrt(px*px+py*py);
    
    } // end of loop over particles
  } // if LHEEventInfo is found


  // first save the vector of candidates
  unsigned int nParticles_=0;
  bool hasStatusFlag=false; // status flags are not always in every sample, need to check first
  std::vector<const reco::Candidate*> cands;
  std::vector<std::vector<reco::GenParticle>::const_iterator> myParticles;
  for( std::vector<reco::GenParticle>::const_iterator 
	 it_gen = genParticleHandle->begin(); 
       it_gen != genParticleHandle->end(); it_gen++ ) 
    {
      reco::GenParticle gen = *it_gen;
      nParticles_++;
       if(nParticles_==1){
       	const reco::GenStatusFlags& cmsswStatus = it_gen->statusFlags();
       	int status=0;
       	for(unsigned int ist=0; ist<cmsswStatus.flags_.size(); ist++)
       	  {
       	    if(cmsswStatus.flags_[ist])
       	      status |= (0x1 << (ist+1)); 
       	  }
       	if(status>0)hasStatusFlag=true;
       	else hasStatusFlag=false;
       }
       // particle needs to have status < 30 or isPrompt if statusFlag information is available
       if( ( applyStatusSelection_ && !applyPromptSelection_ && gen.status()>=30 )
	   || 
	   ( applyPromptSelection_ &&
	     ( 
	      ( hasStatusFlag && !gen.statusFlags().isPrompt() && gen.status()>=30) || 
	       (!hasStatusFlag && gen.status()>=30)
	       )
	     )
	 )continue; // only save beam particle, hard scattering and stable particles
      cands.push_back(&*it_gen);
      myParticles.push_back(it_gen);
    }

  // now loop
  std::vector<const reco::Candidate*>::const_iterator found = cands.begin();
  for(unsigned int genIndex=0; genIndex < MAXNGENPAR_ && genIndex < myParticles.size(); genIndex++){
    
    std::vector<reco::GenParticle>::const_iterator geni = myParticles[genIndex];
    nGenPar_++;

   TLorentzVector p4(geni->px(),geni->py(),geni->pz(),geni->energy());
    new( (*genParP4_)[nGenPar_-1]) TLorentzVector(p4);

    genParQ_.push_back(geni->charge());
    genParId_.push_back(geni->pdgId());
    genParSt_.push_back(geni->status());

    int mompid = -9999;
    if( geni->numberOfMothers() ==1 ) 
      mompid = geni->mother()->pdgId();
    else
      mompid = 10000+geni->numberOfMothers();

    genMomParId_.push_back(mompid);

    genParIndex_.push_back(genIndex);

    int iMo1 = -1;
    int iMo2 = -1;
    int iDa1 = -1;
    int iDa2 = -1;
    int NMo = geni->numberOfMothers();
    int NDa = geni->numberOfDaughters();

    found = find(cands.begin(), cands.end(), geni->mother(0));
    if(found != cands.end()) iMo1 = found - cands.begin() ;

    found = find(cands.begin(), cands.end(), geni->mother(1));
    if(found != cands.end()) iMo2 = found - cands.begin() ;

    found = find(cands.begin(), cands.end(), geni->daughter(0));
    if(found != cands.end()) iDa1 = found - cands.begin() ;

    found = find(cands.begin(), cands.end(), geni->daughter(1));
    if(found != cands.end()) iDa2 = found - cands.begin() ;

    genNMo_.push_back(NMo);
    genNDa_.push_back(NDa);
    genMo1_.push_back(iMo1);
    genMo2_.push_back(iMo2);
    genDa1_.push_back(iDa1);
    genDa2_.push_back(iDa2);
    const reco::GenStatusFlags& cmsswStatus = geni->statusFlags();
    int status=0;
    for(unsigned int ist=0; ist<cmsswStatus.flags_.size(); ist++)
      {
	if(cmsswStatus.flags_[ist])
	  status |= (0x1 << (ist+1)); 
      }
    genStFlag_.push_back(status);
      
  } // end of loop over particles


  //add variables for generator-level study
  Handle<reco::GenMETCollection> metHandle_true;
  if(iEvent.getByToken(genMETToken_true, metHandle_true))
    genMET_true_ = metHandle_true.product()->begin()->pt();
    
  Handle<reco::GenMETCollection> metHandle_calo;
  if(iEvent.getByToken(genMETToken_calo, metHandle_calo))
    genMET_calo_ = metHandle_calo.product()->begin()->pt();
  
  Handle<reco::GenMETCollection> metHandle_caloNonPrompt;
  if(iEvent.getByToken(genMETToken_caloNonPrompt, metHandle_caloNonPrompt))
    genMET_caloNonPrompt_ = metHandle_caloNonPrompt.product()->begin()->pt();

  if(!saveGenJets_)return;

  //ak4genjets
  Handle<reco::GenJetCollection> ak4genJetsHandle;
  if(iEvent.getByToken(ak4genJetsToken,ak4genJetsHandle)){ 
    const reco::GenJetCollection* genJetColl = &(*ak4genJetsHandle);
    reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();   

    for(; gjeti!=genJetColl->end();gjeti++){
	reco::GenJet gjet = *gjeti;
	if(gjet.pt()<=15)continue;
	if(fabs(gjet.eta())>3.0)continue;
	TLorentzVector thisGJet_l4(gjet.px(),gjet.py(),gjet.pz(),gjet.energy());
	new( (*ak4GenJetP4_)[ak4nGenJet_]) TLorentzVector(thisGJet_l4);
	ak4nGenJet_++;
    }
      
  } // end of ak4jet block


  //ak8genjets
  Handle<reco::GenJetCollection> ak8genJetsHandle;
  if(iEvent.getByToken(ak8genJetsToken,ak8genJetsHandle)){ 
    const reco::GenJetCollection* genJetColl = &(*ak8genJetsHandle);
    reco::GenJetCollection::const_iterator gjeti = genJetColl->begin();   

    for(; gjeti!=genJetColl->end();gjeti++){
	reco::GenJet gjet = *gjeti;
	if(gjet.pt()<=100)continue;
	if(fabs(gjet.eta())>3.0)continue;
	TLorentzVector thisGJet_l4(gjet.px(),gjet.py(),gjet.pz(),gjet.energy());
	new( (*ak8GenJetP4_)[ak8nGenJet_]) TLorentzVector(thisGJet_l4);
	ak8nGenJet_++;

	if(!saveGenJetSub_)continue;

	// computing generator-level substructure variables, added by Eiko
	std::vector<const reco::GenParticle*> constituents = gjet.getGenConstituents();
	
	typedef std::vector<fastjet::PseudoJet> VPseudoJet;
	VPseudoJet vjet;
	for(unsigned int ig=0; ig < constituents.size(); ig++){
	  // create vector of PseudoJets
	  const reco::GenParticle *constituent = constituents[ig];
	  if (constituent->pt()<0.01) 
	    continue;
	  vjet.emplace_back(constituent->px(),constituent->py(),constituent->pz(),constituent->energy());
	}
	
	fastjet::ClusterSequenceArea seq(vjet, *jetDefAKT, *areaDef); 
	VPseudoJet alljets = fastjet::sorted_by_pt(seq.inclusive_jets(0.1));

	if (alljets.size()>0){
	  fastjet::PseudoJet *leadingJet = &(alljets[0]);

	  //nsubjettiness
	  
	  ak8GenJettau1_.push_back((*fjtau1)(*leadingJet));
	  ak8GenJettau2_.push_back((*fjtau2)(*leadingJet));
	  ak8GenJettau3_.push_back((*fjtau3)(*leadingJet));


	  fastjet::PseudoJet sdJet = (*softdrop)(*leadingJet);
	  // get and filter constituents of groomed jet
	  VPseudoJet sdconsts = fastjet::sorted_by_pt(sdJet.constituents());

	  ak8GenJetMSD_.push_back(sdJet.m());    
	  ak8GenJetSDSJdR_.push_back(sdJet.structure_of<fastjet::contrib::SoftDrop>().delta_R());
	  ak8GenJetSDSJSymm_.push_back(sdJet.structure_of<fastjet::contrib::SoftDrop>().symmetry());       
	  ak8GenJetSDMassDrop_.push_back(sdJet.structure_of<fastjet::contrib::SoftDrop>().mu());

	} // if find a softdrop jet
	else
	  {

	    ak8GenJetMSD_.push_back(DUMMY);        
	    ak8GenJetSDSJdR_.push_back(DUMMY);     
	    ak8GenJetSDSJSymm_.push_back(DUMMY);   
	    ak8GenJetSDMassDrop_.push_back(DUMMY); 

	  }
    
    } // end of loop over ak8jet  
   
  } // end of ak8jet block

   
}



void  
genInfoTree::SetBranches(){

  AddBranch(&genMET_true_,"genMET_true");
  AddBranch(&genMET_calo_,"genMET_calo");
  AddBranch(&genMET_caloNonPrompt_,"genMET_caloNonPrompt");
  

  AddBranch(&ptHat_, "ptHat");
  AddBranch(&mcWeight_, "mcWeight");

  AddBranch(&HT_, "HT");
  AddBranch(&pdf_, "pdf");
  AddBranch(&originalLHEweight_, "originalLHEweight");
  AddBranch(&pdfscaleSysWeights_, "pdfscaleSysWeights");

  AddBranch(&nGenPar_, "nGenPar");
  AddBranch(&genParP4_, "genParP4");
  AddBranch(&genParQ_,"genParQ");
  AddBranch(&genParId_,"genParId");
  AddBranch(&genParSt_,"genParSt");
  AddBranch(&genMomParId_,"genMomParId");
  AddBranch(&genParIndex_,"genParIndex");

  AddBranch(&genNMo_,"genNMo");
  AddBranch(&genNDa_,"genNDa");
  AddBranch(&genMo1_,"genMo1");
  AddBranch(&genMo2_,"genMo2");
  AddBranch(&genDa1_,"genDa1");
  AddBranch(&genDa2_,"genDa2");
  AddBranch(&genStFlag_,"genStFlag");


  AddBranch(&ak4nGenJet_,  "ak4nGenJet");
  AddBranch(&ak4GenJetP4_, "ak4GenJetP4");

  AddBranch(&ak8nGenJet_,  "ak8nGenJet");
  AddBranch(&ak8GenJetP4_, "ak8GenJetP4");


 /// genjet substructure, added by Eiko
  AddBranch(&ak8GenJetMSD_, "ak8GenJetMSD");       
  AddBranch(&ak8GenJetSDSJdR_, "ak8GenJetSDSJdR");      
  AddBranch(&ak8GenJetSDSJSymm_, "ak8GenJetSDSJSymm");   
  AddBranch(&ak8GenJetSDMassDrop_, "ak8GenJetSDMassDrop"); 
  AddBranch(&ak8GenJettau1_, "ak8GenJettau1");
  AddBranch(&ak8GenJettau2_, "ak8GenJettau2");
  AddBranch(&ak8GenJettau3_, "ak8GenJettau3");
  
  
}


void  
genInfoTree::Clear(){

  ptHat_                = DUMMY;
  mcWeight_             = DUMMY; 
  HT_                   = DUMMY;
  genMET_true_          = DUMMY;
  genMET_calo_          = DUMMY;  
  genMET_caloNonPrompt_ = DUMMY; 

  pdf_.clear();
  originalLHEweight_ = 1;
  pdfscaleSysWeights_.clear();
  nGenPar_ =0;
  genParP4_->Clear();

  genParQ_.clear();
  genParId_.clear();
  genParSt_.clear();
  genMomParId_.clear();
  genParIndex_.clear();
  genNMo_.clear();
  genNDa_.clear();
  genMo1_.clear();
  genMo2_.clear();
  genDa1_.clear();
  genDa2_.clear();
  genStFlag_.clear();

  ak4nGenJet_=0;
  ak4GenJetP4_->Clear();

  ak8nGenJet_=0;
  ak8GenJetP4_->Clear();

 /// genjet substructure, added by Eiko
  ak8GenJetMSD_.clear();       
  ak8GenJetSDSJdR_.clear();      
  ak8GenJetSDSJSymm_.clear();   
  ak8GenJetSDMassDrop_.clear(); 
  ak8GenJettau1_.clear();
  ak8GenJettau2_.clear();
  ak8GenJettau3_.clear();


  
}



