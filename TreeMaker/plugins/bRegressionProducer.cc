#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
// #include "flashgg/DataFormats/interface/Jet.h"
// #include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "RecoJets/JetProducers/interface/PileupJetIdAlgo.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
// #include "flashgg/DataFormats/interface/VertexCandidateMap.h"
// #include "flashgg/DataFormats/interface/Jet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/Math/interface/deltaR.h"


#include <iostream>
#include <string>
#include <vector>

//#include "DNN/TensorFlow/interface/TensorFlow.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

#define debug 0

using namespace std;
using namespace edm;

//namespace flashgg {

class bRegressionProducer : public EDProducer
{

public:
    bRegressionProducer( const ParameterSet & );
    ~bRegressionProducer(){};
    void InitJet();
    void SetNNVectorVar();
    std::vector<float> EvaluateNN();
private:
    void produce( Event &, const EventSetup & ) override;
    //        std::vector<edm::InputTag> inputTagJets_;
    edm::InputTag inputTagJets_;
    EDGetTokenT<View<pat::Jet> > jetToken_;
    edm::EDGetTokenT<double> rhoToken_;
    string bRegressionWeightfileName_;
    double y_mean_,y_std_;

    tensorflow::Session* session;
    tensorflow::GraphDef* graphDef;
    std::vector<float> NNvectorVar_;
    //add vector of mva for eache jet

    //mva variables
    float Jet_pt ;
    float Jet_eta ;
    float rho ;
    float Jet_mt ;
    float Jet_leadTrackPt ;
    float Jet_leptonPtRel ;
    float Jet_leptonDeltaR ;
    float Jet_neHEF ;
    float Jet_neEmEF ;
    float Jet_vtxPt ;
    float Jet_vtxMass ;
    float Jet_vtx3dL ;
    float Jet_vtxNtrk ;
    float Jet_vtx3deL ;
    float Jet_numDaughters_pt03 ;
    float Jet_energyRing_dR0_em_Jet_e ;
    float Jet_energyRing_dR1_em_Jet_e ;
    float Jet_energyRing_dR2_em_Jet_e ;
    float Jet_energyRing_dR3_em_Jet_e ;
    float Jet_energyRing_dR4_em_Jet_e ;
    float Jet_energyRing_dR0_neut_Jet_e ;
    float Jet_energyRing_dR1_neut_Jet_e ;
    float Jet_energyRing_dR2_neut_Jet_e ;
    float Jet_energyRing_dR3_neut_Jet_e ;
    float Jet_energyRing_dR4_neut_Jet_e ;
    float Jet_energyRing_dR0_ch_Jet_e ;
    float Jet_energyRing_dR1_ch_Jet_e ;
    float Jet_energyRing_dR2_ch_Jet_e ;
    float Jet_energyRing_dR3_ch_Jet_e ;
    float Jet_energyRing_dR4_ch_Jet_e ;
    float Jet_energyRing_dR0_mu_Jet_e ;
    float Jet_energyRing_dR1_mu_Jet_e ;
    float Jet_energyRing_dR2_mu_Jet_e ;
    float Jet_energyRing_dR3_mu_Jet_e ;
    float Jet_energyRing_dR4_mu_Jet_e ;
    float Jet_chHEF;//implement from here
    float Jet_chEmEF;
    float Jet_leptonPtRelInv;
    int isEle;
    int isMu;
    int isOther;
    float Jet_mass;
    float Jet_withPtd;

    float CHF_;
    float CEMF_;
    float NHF_ ;
    float NEMF_ ;
    float NumConst_;
    float NumNeutralParticles_;
    float CHM_;
    float MUF_;

    int looseJetID_2016;
    int tightJetID_2017;
    int tightJetID_2018;
};


bRegressionProducer::bRegressionProducer( const ParameterSet &iConfig ) :
    //     inputTagJets_( iConfig.getParameter<std::vector<edm::InputTag> >( "JetTag" )) {
    inputTagJets_( iConfig.getParameter<edm::InputTag>( "JetTag" )) ,
    rhoToken_( consumes<double>(iConfig.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
    //bRegressionWeightfileName_( iConfig.getUntrackedParameter<std::string>("bRegressionWeightfile")),
    y_mean_(iConfig.getUntrackedParameter<double>("y_mean")),
    y_std_(iConfig.getUntrackedParameter<double>("y_std"))
{


    jetToken_= consumes<View<pat::Jet> >(inputTagJets_);

    std::string cmssw_base = getenv("CMSSW_BASE");
    std::string bRegressionWeightfileName_ = cmssw_base+"/src/ExoPieElement/MetaData/data/DNN_models/breg_training_2017.pb";

    //        NNgraph_ = *(new dnn::tf::Graph(bRegressionWeightfileName_.c_str())); //FIXME make this configurable, for variables for breg check this PR https://github.com/cms-analysis/flashgg/pull/968 REMEMBER TO ADD THE LAST CONE!

    graphDef= tensorflow::loadGraphDef(bRegressionWeightfileName_.c_str());
    session = tensorflow::createSession(graphDef);


    Jet_pt = 0.;
    Jet_eta = 0.;
    rho = 0.;
    Jet_mt = 0.;
    Jet_leadTrackPt = 0.;
    Jet_leptonPtRel = 0.;
    Jet_leptonDeltaR = 0.;
    Jet_neHEF = 0.;
    Jet_neEmEF = 0.;
    Jet_vtxPt = 0.;
    Jet_vtxMass = 0.;
    Jet_vtx3dL = 0.;
    Jet_vtxNtrk = 0.;
    Jet_vtx3deL = 0.;
    Jet_numDaughters_pt03 = 0;
    Jet_energyRing_dR0_em_Jet_e = 0.;
    Jet_energyRing_dR1_em_Jet_e = 0.;
    Jet_energyRing_dR2_em_Jet_e = 0.;
    Jet_energyRing_dR3_em_Jet_e = 0.;
    Jet_energyRing_dR4_em_Jet_e = 0.;
    Jet_energyRing_dR0_neut_Jet_e = 0.;
    Jet_energyRing_dR1_neut_Jet_e = 0.;
    Jet_energyRing_dR2_neut_Jet_e = 0.;
    Jet_energyRing_dR3_neut_Jet_e = 0.;
    Jet_energyRing_dR4_neut_Jet_e = 0.;
    Jet_energyRing_dR0_ch_Jet_e = 0.;
    Jet_energyRing_dR1_ch_Jet_e = 0.;
    Jet_energyRing_dR2_ch_Jet_e = 0.;
    Jet_energyRing_dR3_ch_Jet_e = 0.;
    Jet_energyRing_dR4_ch_Jet_e = 0.;
    Jet_energyRing_dR0_mu_Jet_e = 0.;
    Jet_energyRing_dR1_mu_Jet_e = 0.;
    Jet_energyRing_dR2_mu_Jet_e = 0.;
    Jet_energyRing_dR3_mu_Jet_e = 0.;
    Jet_energyRing_dR4_mu_Jet_e = 0.;
    Jet_chHEF = 0.;//implement from here
    Jet_chEmEF = 0.;
    Jet_leptonPtRelInv = 0.;
    isEle = 0.;
    isMu = 0.;
    isOther = 0.;
    Jet_mass = 0.;
    Jet_withPtd = 0.;
    CHF_ = 0;
    CEMF_ = 0;
    NHF_ = 0;
    NEMF_ = 0;
    NumConst_ = 0;
    NumNeutralParticles_ = 0;
    CHM_ = 0;
    MUF_ = 0;

    produces<vector<pat::Jet> > ();
}



void bRegressionProducer::produce( Event &evt, const EventSetup & )
{
    InitJet();
    // input jets
    Handle<View<pat::Jet> > jets;
    evt.getByToken( jetToken_, jets );//just to try get the first one
   unique_ptr<vector<pat::Jet> > jetColl( new vector<pat::Jet> );
    for( unsigned int i = 0 ; i < jets->size() ; i++ ) {


        Ptr<pat::Jet> pjet = jets->ptrAt( i );
        pat::Jet fjet = pat::Jet( *pjet );


        //std::cout << "sv test"  << pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).numberOfSourceCandidatePtrs() << std::endl;
        //            if (fjet.pt()<15. || fabs(fjet.eta())>2.5) continue;


        //variables needed for regression
        Jet_pt = fjet.pt();
        Jet_eta = fjet.eta() ;
        //Jet_leadTrackPt = fjet.userFloat("leadTrackPt");
        edm::Handle<double> rhoHandle;
        evt.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );
        rho = rhoFixedGrd;
        Jet_mt = sqrt(fjet.energy()*fjet.energy()-fjet.pz()*fjet.pz());//seems correct but check again

        // this part we are taking from jetProducer file

         float ptD=0.;
         //int nSecVertices = 0;
         float sumWeight=0;
         float sumPt=0;
         for(const auto & d : pjet->daughterPtrVector()){
             sumWeight+=(d->pt())*(d->pt());
             sumPt+=d->pt();
         }
         ptD = (sumWeight > 0 ? sqrt(sumWeight)/sumPt : 0);

        float leadTrackPt_ = 0, softLepPt = 0, softLepRatio = 0, softLepDr = 0;
        float sumPtDrSq = 0.;
        float sumPtSq = 0.;
        float softLepPtRel = 0.;
        float softLepPtRelInv=0.;


        float cone_boundaries[] = { 0.05, 0.1, 0.2, 0.3, 0.4 }; // hardcoded boundaries: should be made configurable
        size_t ncone_boundaries = sizeof(cone_boundaries)/sizeof(float);
        std::vector<float> chEnergies(ncone_boundaries+1,0.);
        std::vector<float> emEnergies(ncone_boundaries+1,0.);
        std::vector<float> neEnergies(ncone_boundaries+1,0.);
        std::vector<float> muEnergies(ncone_boundaries+1,0.);
        int numDaug03 = 0;

        int softLepPdgId=0;


        for ( unsigned k = 0; k < fjet.numberOfSourceCandidatePtrs(); ++k ) {
            reco::CandidatePtr pfJetConstituent = fjet.sourceCandidatePtr(k);

            const reco::Candidate* kcand = pfJetConstituent.get();
            const pat::PackedCandidate* lPack = dynamic_cast<const pat::PackedCandidate *>( kcand );
            if ( !lPack ) throw cms::Exception( "NoPackedConstituent" ) << " For jet " << i << " failed to get constituent " << k << std::endl;
            float candPt = kcand->pt();
            float candDr   = reco::deltaR(*kcand,fjet);
            sumPtDrSq += candPt*candPt*candDr*candDr;
            sumPtSq += candPt*candPt;

            if( candPt > 0.3 ) { ++numDaug03; }
            if(lPack->charge() != 0 && candPt > leadTrackPt_) leadTrackPt_ = candPt;

            if(abs(lPack->pdgId()) == 11 || abs(lPack->pdgId()) == 13) {
                if(candPt > softLepPt){
                    softLepPt = candPt;
                    softLepRatio = candPt/pjet->pt();
                    if (0){std::cout << "softLepRatio"  << softLepRatio << std::endl;}
                    softLepDr = candDr;
                    softLepPtRel = ( pjet->px()*lPack->px() + pjet->py()*lPack->py() + pjet->pz()*lPack->pz() ) / pjet->p();
                    softLepPtRel = sqrt( lPack->p()*lPack->p() - softLepPtRel*softLepPtRel );

                    softLepPtRelInv = ( pjet->px()*lPack->px() + pjet->py()*lPack->py() + pjet->pz()*lPack->pz() ) / lPack->p();
                    softLepPtRelInv = sqrt( pjet->p()*pjet->p() - softLepPtRelInv*softLepPtRelInv );

                    softLepPdgId = lPack->pdgId();
                }
            }

            int pdgid = abs(lPack->pdgId());
            size_t icone = std::lower_bound(&cone_boundaries[0],&cone_boundaries[ncone_boundaries],candDr) - &cone_boundaries[0];
            float candEnergy = kcand->energy();
            // std::cout << "pdgId " << pdgid << " candDr " << candDr << " icone " << icone << " " << candEnergy << std::endl;
            if( pdgid == 22 || pdgid == 11 ) {
                // std::cout << " fill EM" << std::endl;
                emEnergies[icone] += candEnergy;
            } else if ( pdgid == 13 ) {
                // std::cout << " fill mu" << std::endl;
                muEnergies[icone] += candEnergy;
            } else if ( lPack-> charge() != 0 ) {
                // std::cout << " fill ch" << std::endl;
                chEnergies[icone] += candEnergy;
            } else {
                // std::cout << " fill ne" << std::endl;
                neEnergies[icone] += candEnergy;
            }
        }

        //this max probably not needed, it's just heppy
        Jet_leadTrackPt = leadTrackPt_;
        Jet_leptonPtRel = softLepPtRel;//std::max(float(0.),fjet.userFloat("softLepPtRel"));
        Jet_leptonDeltaR = softLepDr;//std::max(float(0.),fjet.userFloat("softLepDr"));
        Jet_neHEF = fjet.neutralHadronEnergyFraction();
        Jet_neEmEF = fjet.neutralEmEnergyFraction();
        Jet_chHEF = fjet.chargedHadronEnergyFraction();
        Jet_chEmEF = fjet.chargedEmEnergyFraction();
        Jet_leptonPtRelInv = softLepPtRelInv*Jet_pt/fjet.pt();//fjet.userFloat("softLepPtRelInv");
        NHF_  = fjet.neutralHadronEnergyFraction();
        NEMF_ = fjet.neutralEmEnergyFraction();
        CHF_  = fjet.chargedHadronEnergyFraction();
        CEMF_ = fjet.chargedEmEnergyFraction();
        NumConst_ = fjet.chargedMultiplicity()+fjet.neutralMultiplicity();
        NumNeutralParticles_ =fjet.neutralMultiplicity();
        CHM_  = fjet.chargedMultiplicity();
        MUF_  = fjet.muonEnergyFraction();

        int lepPdgID = softLepPdgId;//fjet.userInt("softLepPdgId");
        if (abs(lepPdgID)==13){
            isMu=1;
        }else if (abs(lepPdgID)==11){
            isEle=1;
        }else{
            isOther=1;
        }
        Jet_mass=fjet.mass();
        Jet_withPtd=ptD;//fjet.userFloat("ptD");

        //if(fjet.userFloat("nSecVertices")>0){
//                float vertexX=fjet.userFloat("vtxPosX")-fjet.userFloat("vtxPx");//check if it's correct
//                float vertexY=fjet.userFloat("vtxPosY")-fjet.userFloat("vtxPy");
//                Jet_vtxPt = sqrt(vertexX*vertexX+vertexY*vertexY);
        Jet_vtxPt=0.0;//sqrt(fjet.userFloat("vtxPx")*fjet.userFloat("vtxPx")+fjet.userFloat("vtxPy")*fjet.userFloat("vtxPy"));
        Jet_vtxMass =0.0; //std::max(float(0.),fjet.userFloat("vtxMass"));
        Jet_vtx3dL = 0.0;//std::max(float(0.),fjet.userFloat("vtx3DVal"));
        Jet_vtxNtrk = 0.0;//std::max(float(0.),fjet.userFloat("vtxNTracks"));
        Jet_vtx3deL = 0.0;//std::max(float(0.),fjet.userFloat("vtx3DSig"));
        if (pjet->hasTagInfo("pfSecondaryVertex")){
          float vtxPx=0.0; float vtxPy=0.0; float vtxPz=0.0; float vtxMass=0.0; int vtxNTracks=0; float vtx3DSig=0.0; float vtx3DVal=0.0;
          vtxPx = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).p4().px();
          vtxPy = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).p4().py();
          //vtxPz = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).p4().pz();
          vtxMass = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).p4().mass();
          vtxNTracks = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->secondaryVertex(0).numberOfSourceCandidatePtrs();
          vtx3DSig = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0).significance();
          vtx3DVal = pjet->tagInfoCandSecondaryVertex("pfSecondaryVertex")->flightDistance(0).value();
          Jet_vtxMass=vtxMass;
          Jet_vtx3dL=vtx3DVal;
          Jet_vtxNtrk=vtxNTracks;
          Jet_vtx3deL=vtx3DSig;
          Jet_vtxPt=sqrt(vtxPx*vtxPx+vtxPy*vtxPy);

        }
        //}
        //if (fjet.emEnergies().size()>0){//since in order to save space we save this info only if the candidate has a minimum pt or eta
        Jet_energyRing_dR0_em_Jet_e = emEnergies[0]/fjet.energy();//remember to divide by jet energy
        Jet_energyRing_dR1_em_Jet_e = emEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_em_Jet_e = emEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_em_Jet_e = emEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_em_Jet_e = emEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_neut_Jet_e = neEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_neut_Jet_e = neEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_neut_Jet_e = neEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_neut_Jet_e = neEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_neut_Jet_e = neEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_ch_Jet_e = chEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_ch_Jet_e = chEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_ch_Jet_e = chEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_ch_Jet_e = chEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_ch_Jet_e = chEnergies[4]/fjet.energy();
        Jet_energyRing_dR0_mu_Jet_e = muEnergies[0]/fjet.energy();
        Jet_energyRing_dR1_mu_Jet_e = muEnergies[1]/fjet.energy();
        Jet_energyRing_dR2_mu_Jet_e = muEnergies[2]/fjet.energy();
        Jet_energyRing_dR3_mu_Jet_e = muEnergies[3]/fjet.energy();
        Jet_energyRing_dR4_mu_Jet_e = muEnergies[4]/fjet.energy();
        //}
        Jet_numDaughters_pt03 = numDaug03;//fjet.userInt("numDaug03");

        std::vector<float> bRegNN(3,-999);


        if(debug){
            cout<<"Jet_pt :"<<Jet_pt <<endl;
            cout<<"Jet_eta :"<<Jet_eta <<endl;
            cout<<"rho :"<<rho <<endl;
            cout<<"Jet_mt :"<<Jet_mt <<endl;
            cout<<"Jet_leadTrackPt :"<<Jet_leadTrackPt <<endl;
            cout<<"Jet_leptonPtRel :"<<Jet_leptonPtRel <<endl;
            cout<<"Jet_leptonDeltaR :"<<Jet_leptonDeltaR <<endl;
            cout<<"Jet_neHEF :"<<Jet_neHEF <<endl;
            cout<<"Jet_neEmEF :"<<Jet_neEmEF <<endl;
            cout<<"Jet_vtxPt :"<<Jet_vtxPt <<endl;
            cout<<"Jet_vtxMass :"<<Jet_vtxMass <<endl;
            cout<<"Jet_vtx3dL :"<<Jet_vtx3dL <<endl;
            cout<<"Jet_vtxNtrk :"<<Jet_vtxNtrk <<endl;
            cout<<"Jet_vtx3deL :"<<Jet_vtx3deL <<endl;
            cout<<"Jet_energyRing_dR0_em_Jet_e :"<<Jet_energyRing_dR0_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_em_Jet_e :"<<Jet_energyRing_dR1_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_em_Jet_e :"<<Jet_energyRing_dR2_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_em_Jet_e :"<<Jet_energyRing_dR3_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_em_Jet_e :"<<Jet_energyRing_dR4_em_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_neut_Jet_e :"<<Jet_energyRing_dR0_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_neut_Jet_e :"<<Jet_energyRing_dR1_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_neut_Jet_e :"<<Jet_energyRing_dR2_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_neut_Jet_e :"<<Jet_energyRing_dR3_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_neut_Jet_e :"<<Jet_energyRing_dR4_neut_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_ch_Jet_e :"<<Jet_energyRing_dR0_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_ch_Jet_e :"<<Jet_energyRing_dR1_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_ch_Jet_e :"<<Jet_energyRing_dR2_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_ch_Jet_e :"<<Jet_energyRing_dR3_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_ch_Jet_e :"<<Jet_energyRing_dR4_ch_Jet_e <<endl;
            cout<<"Jet_energyRing_dR0_mu_Jet_e :"<<Jet_energyRing_dR0_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR1_mu_Jet_e :"<<Jet_energyRing_dR1_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR2_mu_Jet_e :"<<Jet_energyRing_dR2_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR3_mu_Jet_e :"<<Jet_energyRing_dR3_mu_Jet_e <<endl;
            cout<<"Jet_energyRing_dR4_mu_Jet_e :"<<Jet_energyRing_dR4_mu_Jet_e <<endl;
            cout<<"Jet_numDaughters_pt03 :"<<Jet_numDaughters_pt03 <<endl;



        }

        SetNNVectorVar();
        bRegNN = EvaluateNN();
        NNvectorVar_.clear();

        //FIXME read through file, config is here /afs/cern.ch/user/n/nchernya/public/100M_2018-03-01_job23_rawJetsJECtarget/config.json
        //            float y_mean= 1.0454729795455933;
        //            float y_std = 0.3162831664085388;

        fjet.addUserFloat("bRegNNCorr", bRegNN[0]*y_std_+y_mean_);
        fjet.addUserFloat("bRegNNResolution",0.5*(bRegNN[2]-bRegNN[1])*y_std_);

        //std::cout << "bRegNNCorr" << bRegNN[0]*y_std_+y_mean_  << std::endl;
        //std::cout << "checking userfloat" << fjet.userFloat("bRegNNCorr")  << std::endl;

        if  (abs(Jet_eta)<=2.7) {
            looseJetID_2016 = (NHF_<0.99 && NEMF_<0.99 && NumConst_>1) && ((abs(Jet_eta)<=2.4 && CHF_>0 && CHM_>0 && CEMF_<0.99) || abs(Jet_eta)>2.4) && abs(Jet_eta)<=2.7;
            tightJetID_2017 = (NHF_<0.90 && NEMF_<0.90 && NumConst_>1) && ((abs(Jet_eta)<=2.4 && CHF_>0 && CHM_>0) || abs(Jet_eta)>2.4) && abs(Jet_eta)<=2.7;
            if  (abs(Jet_eta)<=2.6) {
                tightJetID_2018=(abs(eta)<=2.6 && CEMF_<0.8 && CHM_>0 && CHF_>0 && NumConst_>1 && NEMF_<0.9 && MUF_ <0.8 && NHF_ < 0.9 );
            }
            if  (abs(Jet_eta)>2.6 && abs(Jet_eta)<= 2.7){
                tightJetID_2018=( abs(eta)>2.6 && abs(eta)<=2.7 && CEMF_<0.8 && CHM_>0 && NEMF_<0.99 && MUF_ <0.8 && NHF_ < 0.9 );
            }
        }
        else if  (abs(Jet_eta)>2.7 && abs(Jet_eta)<= 3.0){
            looseJetID_2016 = (NHF_<0.98 && NEMF_>0.01 && NumNeutralParticles_>2 && abs(Jet_eta)>2.7 && abs(Jet_eta)<=3.0 );
            tightJetID_2017 = (NEMF_<0.99 && NEMF_>0.02 && NumNeutralParticles_>2 && abs(Jet_eta)>2.7 && abs(Jet_eta)<=3.0 );
            tightJetID_2018 = (NEMF_>0.02 && NEMF_<0.99 && NumNeutralParticle_>2 && abs(eta)>2.7 && abs(eta)<=3.0 )
        }
        else if (abs(Jet_eta)>3.0) {
            looseJetID_2016 = (NEMF_<0.90 && NumNeutralParticles_>10 && abs(Jet_eta)>3.0 );
            tightJetID_2017 = (NEMF_<0.90 && NHF_ >0.02 && NumNeutralParticles_>10 && abs(Jet_eta)>3.0 );
            tightJetID_2018 = (NEMF<0.90 && NHF>0.2 && NumNeutralParticle>10 && abs(eta)>3.0 )
        }
        fjet.addUserInt("looseJetID_2016", looseJetID_2016);
        fjet.addUserInt("tightJetID_2017", tightJetID_2017);
        fjet.addUserFloat("NHF_",NHF_);
        fjet.addUserFloat("NEMF_",NEMF_);
        fjet.addUserFloat("CHF_",CHF_);
        fjet.addUserFloat("CEMF_",CEMF_);
        fjet.addUserFloat("CHM_",CHM_);
        fjet.addUserFloat("MUF_",MUF_);
        fjet.addUserFloat("NumConst_",NumConst_);
        fjet.addUserFloat("NumNeutralParticles_",NumNeutralParticles_);


        jetColl->push_back( fjet );


    }
    evt.put( std::move( jetColl ) );
}

void bRegressionProducer::InitJet(){
    Jet_pt = 0.;
    Jet_eta = 0.;
    rho = 0.;
    Jet_mt = 0.;
    Jet_leadTrackPt = 0.;
    Jet_leptonPtRel = 0.;
    Jet_leptonDeltaR = 0.;
    Jet_neHEF = 0.;
    Jet_neEmEF = 0.;
    Jet_vtxPt = 0.;
    Jet_vtxMass = 0.;
    Jet_vtx3dL = 0.;
    Jet_vtxNtrk = 0.;
    Jet_vtx3deL = 0.;
    Jet_numDaughters_pt03 = 0;
    Jet_energyRing_dR0_em_Jet_e = 0.;
    Jet_energyRing_dR1_em_Jet_e = 0.;
    Jet_energyRing_dR2_em_Jet_e = 0.;
    Jet_energyRing_dR3_em_Jet_e = 0.;
    Jet_energyRing_dR4_em_Jet_e = 0.;
    Jet_energyRing_dR0_neut_Jet_e = 0.;
    Jet_energyRing_dR1_neut_Jet_e = 0.;
    Jet_energyRing_dR2_neut_Jet_e = 0.;
    Jet_energyRing_dR3_neut_Jet_e = 0.;
    Jet_energyRing_dR4_neut_Jet_e = 0.;
    Jet_energyRing_dR0_ch_Jet_e = 0.;
    Jet_energyRing_dR1_ch_Jet_e = 0.;
    Jet_energyRing_dR2_ch_Jet_e = 0.;
    Jet_energyRing_dR3_ch_Jet_e = 0.;
    Jet_energyRing_dR4_ch_Jet_e = 0.;
    Jet_energyRing_dR0_mu_Jet_e = 0.;
    Jet_energyRing_dR1_mu_Jet_e = 0.;
    Jet_energyRing_dR2_mu_Jet_e = 0.;
    Jet_energyRing_dR3_mu_Jet_e = 0.;
    Jet_energyRing_dR4_mu_Jet_e = 0.;
    Jet_chHEF = 0.;//implement from here
    Jet_chEmEF = 0.;
    Jet_leptonPtRelInv = 0.;
    isEle = 0.;
    isMu = 0.;
    isOther = 0.;
    Jet_mass = 0.;
    Jet_withPtd = 0.;


}//end InitJet

void bRegressionProducer::SetNNVectorVar(){

    NNvectorVar_.push_back(Jet_pt) ;//0
    NNvectorVar_.push_back(Jet_eta) ;
    NNvectorVar_.push_back(rho) ;
    NNvectorVar_.push_back(Jet_mt) ;
    NNvectorVar_.push_back(Jet_leadTrackPt) ;
    NNvectorVar_.push_back(Jet_leptonPtRel) ;//5
    NNvectorVar_.push_back(Jet_leptonDeltaR) ;
    NNvectorVar_.push_back(Jet_neHEF) ;
    NNvectorVar_.push_back(Jet_neEmEF) ;
    NNvectorVar_.push_back(Jet_vtxPt) ;
    NNvectorVar_.push_back(Jet_vtxMass) ;//10
    NNvectorVar_.push_back(Jet_vtx3dL) ;
    NNvectorVar_.push_back(Jet_vtxNtrk) ;
    NNvectorVar_.push_back(Jet_vtx3deL) ;
    NNvectorVar_.push_back(Jet_numDaughters_pt03) ;//this variable has changed order, in bdt it was last, check why
    NNvectorVar_.push_back(Jet_energyRing_dR0_em_Jet_e) ;//15
    NNvectorVar_.push_back(Jet_energyRing_dR1_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_em_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_neut_Jet_e) ;//20
    NNvectorVar_.push_back(Jet_energyRing_dR1_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_neut_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_ch_Jet_e) ;//25
    NNvectorVar_.push_back(Jet_energyRing_dR1_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_ch_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR0_mu_Jet_e) ;//30
    NNvectorVar_.push_back(Jet_energyRing_dR1_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR2_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR3_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_energyRing_dR4_mu_Jet_e) ;
    NNvectorVar_.push_back(Jet_chHEF);//35
    NNvectorVar_.push_back(Jet_chEmEF);
    NNvectorVar_.push_back(Jet_leptonPtRelInv);
    NNvectorVar_.push_back(isEle);
    NNvectorVar_.push_back(isMu);
    NNvectorVar_.push_back(isOther);//40
    NNvectorVar_.push_back(Jet_mass);
    NNvectorVar_.push_back(Jet_withPtd);


}

std::vector<float> bRegressionProducer::EvaluateNN(){
    tensorflow::Tensor input(tensorflow::DT_FLOAT, {1,43});
    for (unsigned int i = 0; i < NNvectorVar_.size(); i++){
        //            std::cout<<"i:"<<i<<" x:"<<NNvectorVar_[i]<<std::endl;
        input.matrix<float>()(0,i) =  float(NNvectorVar_[i]);
    }
    std::vector<tensorflow::Tensor> outputs;
    tensorflow::run(session, { { "ffwd_inp:0",input } }, { "ffwd_out/BiasAdd:0" }, &outputs);

    std::vector<float> correction(3);//3 outputs, first value is mean and then other 2 quantiles
    correction[0] = outputs[0].matrix<float>()(0, 0);
    correction[1] = outputs[0].matrix<float>()(0, 1);
    correction[2] = outputs[0].matrix<float>()(0, 2);

    //        std::cout<<correction[0]<<" "<<correction[1]<<" "<<correction[2]<<std::endl;

    return correction;

}//end EvaluateNN

//}



//typedef pat::bRegressionProducer flashggbRegressionProducer;
DEFINE_FWK_MODULE( bRegressionProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
