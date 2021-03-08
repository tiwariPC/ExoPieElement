#include "ExoPieElement/TreeMaker/interface/genInfoTree.h"
#include <bitset>
#include "DataFormats/HepMCCandidate/interface/GenStatusFlags.h"
#include "TLorentzVector.h"

genInfoTree::genInfoTree(std::string name, TTree* tree, const edm::ParameterSet& iConfig):
  baseTree(name,tree),
  MAXNGENPAR_(iConfig.getParameter<unsigned int>("maxNumGenPar")),
  applyStatusSelection_(iConfig.getParameter<bool>("applyStatusSelection")),
  applyPromptSelection_(iConfig.getParameter<bool>("applyPromptSelection")),
  saveLHEWeights_(iConfig.getParameter<bool>("saveLHEWeights")),
  saveGenJets_(iConfig.getParameter<bool>("saveGenJets"))
{
  //genParP4_ =   new TClonesArray("TLorentzVector");
  //ak4GenJetP4_ =   new TClonesArray("TLorentzVector");
  //ak8GenJetP4_ =   new TClonesArray("TLorentzVector");

  SetBranches();
}


genInfoTree::~genInfoTree()
{
  //delete genParP4_;
  //delete ak4GenJetP4_;
  //delete ak8GenJetP4_;

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
      for (unsigned int i=0; i< evt->weights().size(); i++)
      {
        float tempMCWeight = genEventInfoHandle->weight() > 0? 1: -1;
        float sysLHEweight = tempMCWeight* evt->weights()[i].wgt/evt->originalXWGTUP();
        pdfscaleSysWeights_.push_back( sysLHEweight );
        pdfscaleSysWgtID_.push_back(evt->weights()[i].id);
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

    //TLorentzVector p4(geni->px(),geni->py(),geni->pz(),geni->energy());
    //new( (*genParP4_)[nGenPar_-1]) TLorentzVector(p4);

    genParPx_.push_back(geni->px());
    genParPy_.push_back(geni->py());
    genParPz_.push_back(geni->pz());
    genParE_.push_back(geni->energy());

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
  Handle<pat::PackedGenParticleCollection> genMETParticleHandle;
  if(not iEvent.getByToken(genMETToken_true, genMETParticleHandle))
  {
      std::cout<<
         "GenAnalyzer: Generator Level Information not found\n"
         <<std::endl;
  }
    TLorentzVector vV;
    std::vector<const pat::PackedCandidate*> candis;
    std::vector<std::vector<pat::PackedGenParticle>::const_iterator> myMETParticles;
    for( std::vector<pat::PackedGenParticle>::const_iterator it_gen = genMETParticleHandle->begin(); it_gen != genMETParticleHandle->end(); it_gen++ )    {
        pat::PackedGenParticle gen = *it_gen;
        //  std::cout<<" px = "<<gen.px()<<std::endl;
      if (abs(gen.pdgId())==18){
          TLorentzVector tmp_;
          tmp_.SetPxPyPzE(gen.px(), gen.py(), gen.pz(), gen.energy());
          vV += tmp_;
          //vV.SetPxPyPzE(gen.px(), gen.py(), gen.pz(), gen.energy());
          //std::cout<<" inside dm  "<<gen.pt()<<" " <<gen.status()<<std::endl;
        }
    }
  genMET_true_ = vV.Pt();

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

	//TLorentzVector thisGJet_l4(gjet.px(),gjet.py(),gjet.pz(),gjet.energy());
	//new( (*ak4GenJetP4_)[ak4nGenJet_]) TLorentzVector(thisGJet_l4);
	ak4GenJetPx_.push_back(gjet.px());
	ak4GenJetPy_.push_back(gjet.py());
	ak4GenJetPz_.push_back(gjet.pz());
	ak4GenJetE_.push_back(gjet.energy());

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
	//TLorentzVector thisGJet_l4(gjet.px(),gjet.py(),gjet.pz(),gjet.energy());
	//new( (*ak8GenJetP4_)[ak8nGenJet_]) TLorentzVector(thisGJet_l4);
	ak8GenJetPx_.push_back(gjet.px());
	ak8GenJetPy_.push_back(gjet.py());
	ak8GenJetPz_.push_back(gjet.pz());
	ak8GenJetE_.push_back(gjet.energy());
	ak8nGenJet_++;
    }
  } // end of ak8jet block
}

bool gen_extra=false;
void
genInfoTree::SetBranches(){

  AddBranch(&nGenPar_, "nGenPar");
  //AddBranch(&genParP4_, "genParP4");
  AddBranch(&genParPx_, "genParPx");
  AddBranch(&genParPy_, "genParPy");
  AddBranch(&genParPz_, "genParPz");
  AddBranch(&genParE_, "genParE");
  AddBranch(&genParId_,"genParId");
  AddBranch(&genParSt_,"genParSt");
  AddBranch(&genMomParId_,"genMomParId");
  AddBranch(&mcWeight_, "mcWeight");
  AddBranch(&originalLHEweight_, "originalLHEweight");
  AddBranch(&pdfscaleSysWgtID_, "pdfscaleSysWgtID_");
  AddBranch(&pdfscaleSysWeights_, "pdfscaleSysWeights");

  if (gen_extra){
    AddBranch(&genParIndex_,"genParIndex");
    AddBranch(&genParQ_,"genParQ");
    AddBranch(&genMET_true_,"genMET_true");
    AddBranch(&genMET_calo_,"genMET_calo");
    AddBranch(&genMET_caloNonPrompt_,"genMET_caloNonPrompt");
    AddBranch(&ptHat_, "ptHat");
    AddBranch(&HT_, "HT");
    AddBranch(&pdf_, "pdf");


    AddBranch(&genNMo_,"genNMo");
    AddBranch(&genNDa_,"genNDa");
    AddBranch(&genMo1_,"genMo1");
    AddBranch(&genMo2_,"genMo2");
    AddBranch(&genDa1_,"genDa1");
    AddBranch(&genDa2_,"genDa2");
    AddBranch(&genStFlag_,"genStFlag");

    AddBranch(&ak4nGenJet_,  "ak4nGenJet");
    //AddBranch(&ak4GenJetP4_, "ak4GenJetP4");
    AddBranch(&ak4GenJetPx_, "ak4GenJetPx");
    AddBranch(&ak4GenJetPy_, "ak4GenJetPy");
    AddBranch(&ak4GenJetPz_, "ak4GenJetPz");
    AddBranch(&ak4GenJetE_, "ak4GenJetE");

    AddBranch(&ak8nGenJet_,  "ak8nGenJet");
    //AddBranch(&ak8GenJetP4_, "ak8GenJetP4");
    AddBranch(&ak8GenJetPx_, "ak8GenJetPx");
    AddBranch(&ak8GenJetPy_, "ak8GenJetPy");
    AddBranch(&ak8GenJetPz_, "ak8GenJetPz");
    AddBranch(&ak8GenJetE_, "ak8GenJetE");
  }
}


void
genInfoTree::Clear(){

  ptHat_                = -9999.0;
  mcWeight_             = -9999.0;
  HT_                   = -9999.0;
  genMET_true_          = -9999.0;
  genMET_calo_          = -9999.0;
  genMET_caloNonPrompt_ = -9999.0;

  pdf_.clear();
  originalLHEweight_ = 1;
  pdfscaleSysWeights_.clear();
  pdfscaleSysWgtID_.clear();
  nGenPar_ =0;
  //genParP4_->Clear();
  genParPx_.clear();
  genParPy_.clear();
  genParPz_.clear();
  genParE_.clear();

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
  //ak4GenJetP4_->Clear();
  ak4GenJetPx_.clear();
  ak4GenJetPy_.clear();
  ak4GenJetPz_.clear();
  ak4GenJetE_.clear();

  ak8nGenJet_=0;
  //ak8GenJetP4_->Clear();
  ak8GenJetPx_.clear();
  ak8GenJetPy_.clear();
  ak8GenJetPz_.clear();
  ak8GenJetE_.clear();



}
