#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <iostream>
#include <string>
#include <vector>

#define debug 0

using namespace std;
using namespace edm;

//namespace flashgg {

class FatJetIDProducer : public EDProducer
{

public:
    FatJetIDProducer( const ParameterSet & );
    ~FatJetIDProducer(){};
    void InitJet();
private:
    void produce( Event &, const EventSetup & ) override;
    edm::InputTag inputTagJets_;
    EDGetTokenT<View<pat::Jet> > jetToken_;
    float Jet_pt ;
    float Jet_eta ;
    float CHF_;
    float CEMF_;
    float NHF_ ;
    float NEMF_ ;
    float NumConst_;
    float NumNeutralParticles_;
    float CHM_;

    int looseJetID_2016;
    int tightJetID_2017;
};


FatJetIDProducer::FatJetIDProducer( const ParameterSet &iConfig ) :
    inputTagJets_( iConfig.getParameter<edm::InputTag>( "JetTag" ))
{
    jetToken_= consumes<View<pat::Jet> >(inputTagJets_);
    Jet_pt = 0.;
    Jet_eta = 0.;
    CHF_ = 0;
    CEMF_ = 0;
    NHF_ = 0;
    NEMF_ = 0;
    NumConst_ = 0;
    NumNeutralParticles_ = 0;
    CHM_ = 0;
    produces<vector<pat::Jet> > ();
}



void FatJetIDProducer::produce( Event &evt, const EventSetup & )
{
    InitJet();
    // input jets
    Handle<View<pat::Jet> > jets;
    evt.getByToken( jetToken_, jets );//just to try get the first one
   unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> );
    for( unsigned int i = 0 ; i < jets->size() ; i++ ) {
        Ptr<pat::Jet> pjet = jets->ptrAt( i );
        pat::Jet fjet = pat::Jet( *pjet );
        //variables needed for regression
        Jet_pt = fjet.pt();
        Jet_eta = fjet.eta() ;

        NHF_  = fjet.neutralHadronEnergyFraction();
        NEMF_ = fjet.neutralEmEnergyFraction();
        CHF_  = fjet.chargedHadronEnergyFraction();
        CEMF_ = fjet.chargedEmEnergyFraction();
        NumConst_ = fjet.chargedMultiplicity()+fjet.neutralMultiplicity();
        NumNeutralParticles_ =fjet.neutralMultiplicity();
        CHM_  = fjet.chargedMultiplicity();


        if(debug){
            cout<<"Jet_pt :"<<Jet_pt <<endl;
            cout<<"Jet_eta :"<<Jet_eta <<endl;
        }
        if  (abs(Jet_eta)<=2.7) {
            looseJetID_2016 = (NHF_<0.99 && NEMF_<0.99 && NumConst_>1) && ((abs(Jet_eta)<=2.4 && CHF_>0 && CHM_>0 && CEMF_<0.99) || abs(Jet_eta)>2.4) && abs(Jet_eta)<=2.7;
            tightJetID_2017 = (NHF_<0.90 && NEMF_<0.90 && NumConst_>1) && ((abs(Jet_eta)<=2.4 && CHF_>0 && CHM_>0) || abs(Jet_eta)>2.4) && abs(Jet_eta)<=2.7;
        }
        else if  (abs(Jet_eta)>2.7 && abs(Jet_eta)<= 3.0){
            looseJetID_2016 = (NHF_<0.98 && NEMF_>0.01 && NumNeutralParticles_>2 && abs(Jet_eta)>2.7 && abs(Jet_eta)<=3.0 );
            tightJetID_2017 = (NEMF_<0.99 && NEMF_>0.02 && NumNeutralParticles_>2 && abs(Jet_eta)>2.7 && abs(Jet_eta)<=3.0 );
        }
        else if (abs(Jet_eta)<=3.0) {
            looseJetID_2016 = (NEMF_<0.90 && NumNeutralParticles_>10 && abs(Jet_eta)>3.0 );
            tightJetID_2017 = (NEMF_<0.90 && NHF_ >0.02 && NumNeutralParticles_>10 && abs(Jet_eta)>3.0 );
        }
        fjet.addUserInt("looseJetID_2016", looseJetID_2016);
        fjet.addUserInt("tightJetID_2017", tightJetID_2017);

        jetColl->push_back( fjet );
    }
    evt.put( std::move( jetColl ) );
}

void FatJetIDProducer::InitJet(){
    Jet_pt = 0.;
    Jet_eta = 0.;

}
//typedef pat::FatJetIDProducer flashggFatJetIDProducer;
DEFINE_FWK_MODULE( FatJetIDProducer );
