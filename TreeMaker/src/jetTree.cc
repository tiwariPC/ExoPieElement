
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
const math::XYZPoint & position(const reco::Vertex & sv) {
        return sv.position();
}
const math::XYZPoint & position(const reco::VertexCompositePtrCandidate & sv) {
        return sv.vertex();
}

jetTree::jetTree(std::string desc, TTree* tree, const edm::ParameterSet& iConfig) :


        baseTree(desc, tree),
        isTHINJet_(false),
        isFATJet_(false),
        isAK4PuppiJet_(false),
        isAK8PuppiJet_(false),
        isCA15PuppiJet_(false),
        useJECText_(iConfig.getParameter<bool>("useJECText")),
        runOn2018_(iConfig.getParameter<bool>("runOn2018")),
        runOn2017_(iConfig.getParameter<bool>("runOn2017")),
        runOn2016_(iConfig.getParameter<bool>("runOn2016")),
        svTagInfosCstr_(iConfig.getParameter<std::string>("svTagInfosPY")),
        jecUncPayLoadName_(iConfig.getParameter<std::string>(Form("%sjecUncPayLoad",desc.data()))),
        jecNames_(iConfig.getParameter<std::vector<std::string> >(Form("%sjecNames",desc.data()) )),
        jecUncName_(iConfig.getParameter<std::string>(Form("%sjecUncName",desc.data())) ),
        jet2017ID_()
{

        if (desc.find("THIN")!=std::string::npos)
                isTHINJet_=true;
        if (desc.find("FAT")!=std::string::npos)
                isFATJet_=true;
        if (desc.find("AK4Puppi")!=std::string::npos)
                isAK4PuppiJet_=true;
        if (desc.find("AK8Puppi")!=std::string::npos)
                isAK8PuppiJet_=true;
        if (desc.find("CA15Puppi")!=std::string::npos)
                isCA15PuppiJet_=true;

        std::cout << " inside jet tree "<< desc << std::endl;


        //genjetP4_    = new TClonesArray("TLorentzVector");
        //jetP4_       = new TClonesArray("TLorentzVector");
        //unCorrJetP4_ = new TClonesArray("TLorentzVector");
        //jetCHSP4_    = new TClonesArray("TLorentzVector");
        jetSDRawP4_  = new TClonesArray("TLorentzVector");

        SetBranches();


        if(isFATJet_)
        {
                prunedMassJecNames_ = iConfig.getParameter<std::vector<std::string> >(Form("%sprunedMassJecNames",desc.data()));

                if(useJECText_) {

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

        std::cout<<" before ca15 "<<std::endl;
        if(isCA15PuppiJet_) {
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

                /* ECF: Ends here */



                std::string cmssw_base = getenv("CMSSW_BASE");
                std::string fweight = cmssw_base+"/src/ExoPieElement/TreeMaker/data/BoostedSVDoubleCA15_withSubjet_v4.weights.xml";
                mJetBoostedBtaggingMVACalc.initialize("BDT",fweight);
        } // if it's CA15Puppijet

        std::cout<<" after ca15 dbt "<<std::endl;
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

        // JEC Unceratinty sources
        const int nsrc = 11;
        const char* srcnames_2016[nsrc] = {"Absolute", "Absolute_2016", "BBEC1", "BBEC1_2016", "EC2", "EC2_2016", "FlavorQCD", "HF", "HF_2016", "RelativeBal", "RelativeSample_2016"};
        const char* srcnames_2017[nsrc] = {"Absolute", "Absolute_2017", "BBEC1", "BBEC1_2017", "EC2", "EC2_2017", "FlavorQCD", "HF", "HF_2017", "RelativeBal", "RelativeSample_2017"};
        const char* srcnames_2018[nsrc] = {"Absolute", "Absolute_2018", "BBEC1", "BBEC1_2018", "EC2", "EC2_2018", "FlavorQCD", "HF", "HF_2018", "RelativeBal", "RelativeSample_2018"};

        std::string cmssw_base = getenv("CMSSW_BASE");
        //std::vector<JetCorrectionUncertainty*> vsrc(nsrc);
        if (runOn2016_) {
                std::string JecSourceUncFile = cmssw_base+"/src/ExoPieElement/TreeMaker/data/RegroupedV2_Summer16_07Aug2017_V11_MC_UncertaintySources_AK4PFchs.txt";
                total = new JetCorrectionUncertainty(*(new JetCorrectorParameters(JecSourceUncFile, "Total")));
                for (int isrc = 0; isrc < nsrc; isrc++) {
                        const char *name = srcnames_2016[isrc];
                        p = new JetCorrectorParameters(JecSourceUncFile, name);
                        unc = new JetCorrectionUncertainty(*p);
                        vsrc.push_back(unc);//vsrc[isrc] = unc;
                };
        };

        if (runOn2017_) {
                std::string JecSourceUncFile = cmssw_base+"/src/ExoPieElement/TreeMaker/data/RegroupedV2_Fall17_17Nov2017_V32_MC_UncertaintySources_AK4PFchs.txt";
                total = new JetCorrectionUncertainty(*(new JetCorrectorParameters(JecSourceUncFile, "Total")));
                for (int isrc = 0; isrc < nsrc; isrc++) {
                        const char *name = srcnames_2017[isrc];
                        p = new JetCorrectorParameters(JecSourceUncFile, name);
                        unc = new JetCorrectionUncertainty(*p);
                        vsrc.push_back(unc);//vsrc[isrc] = unc;
                };
        };
        if (runOn2018_) {
                std::string JecSourceUncFile = cmssw_base+"/src/ExoPieElement/TreeMaker/data/RegroupedV2_Autumn18_V19_MC_UncertaintySources_AK4PFchs.txt";
                total = new JetCorrectionUncertainty(*(new JetCorrectorParameters(JecSourceUncFile, "Total")));
                for (int isrc = 0; isrc < nsrc; isrc++) {
                        const char *name = srcnames_2018[isrc];
                        p = new JetCorrectorParameters(JecSourceUncFile, name);
                        unc = new JetCorrectionUncertainty(*p);
                        vsrc.push_back(unc);//vsrc[isrc] = unc;
                };
        };
        std::cout<<" after jec "<<std::endl;
}


jetTree::~jetTree(){

        //delete genjetP4_;
        //delete jetP4_;
        //delete unCorrJetP4_;
        //delete jetCHSP4_;
        delete jetSDRawP4_;

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




        // Get the primary vertex collection
        edm::Handle<reco::VertexCollection>  h_pv;
        if(not iEvent.getByToken(vertexToken,h_pv))
        {
                std::cout<<"FATAL EXCEPTION: in beginging "<<"Following Not Found: "
                         <<"vertexToken"<<std::endl;
                exit(0);
        }

        if (h_pv->empty()) return; // skip the event if no PV found

        jetNPV_=  h_pv->size();


        // Get the Handle for Jet Collection
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
        if(!useJECText_) {
                edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
                iSetup.get<JetCorrectionsRecord>().get(jecUncPayLoadName_.data(),JetCorParColl);
                JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
                jecUnc_ = new JetCorrectionUncertainty(JetCorPar);
        }



        // Loop over jet collection based on the jet type flag
        pat::JetCollection jets(*(JetHandle.product()));
        std::sort(jets.begin(),jets.end(),PtGreater());
        std::vector<pat::Jet>::const_iterator jet =jets.begin();


        for(; jet!=jets.end(); jet++) {

                if(jet->pt() < 30. || (!(jet->isPFJet()))) continue;
                if((isFATJet_ || isAK8PuppiJet_ || isCA15PuppiJet_) && jet->pt() < 200.) continue;

                nJet_++;
                //Stuff common for all jets.


                if (jet->genJet()) {

                        double DR =0;
                        double a = ( jet->p4().eta() - jet->genJet()->p4().eta()) * ( jet->p4().eta() - jet->genJet()->p4().eta());
                        double b = ( jet->p4().phi() - jet->genJet()->p4().phi()) * ( jet->p4().phi() - jet->genJet()->p4().phi());
                        DR = sqrt(a+b);

                        /*
                           new( (*genjetP4_)[nJet_-1]) TLorentzVector(
                                 jet->genJet()->p4().px(),
                                 jet->genJet()->p4().py(),
                                 jet->genJet()->p4().pz(),
                                 jet->genJet()->p4().energy()
                                 );
                         */
                        genjetpx_.push_back(jet->genJet()->p4().px());
                        genjetpy_.push_back(jet->genJet()->p4().py());
                        genjetpz_.push_back(jet->genJet()->p4().pz());
                        genjetE_.push_back(jet->genJet()->p4().energy());

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
                        genjetpx_.push_back(DUMMY);
                        genjetpy_.push_back(DUMMY);
                        genjetpz_.push_back(DUMMY);
                        genjetE_.push_back(DUMMY);
                        genjetEM_.push_back(DUMMY);
                        genjetHAD_.push_back(DUMMY);
                        genjetINV_.push_back(DUMMY);
                        genjetAUX_.push_back(DUMMY);
                        matchedDR_.push_back(DUMMY);
                        //new( (*genjetP4_)[nJet_-1]) TLorentzVector(DUMMY,DUMMY,DUMMY,DUMMY);
                }

                jetRawFactor_.push_back(jet->jecFactor("Uncorrected"));
                reco::Candidate::LorentzVector uncorrJet;
                uncorrJet = jet->correctedP4(0);
                /*
                   new( (*unCorrJetP4_)[nJet_-1]) TLorentzVector(
                          uncorrJet.px(),
                          uncorrJet.py(),
                          uncorrJet.pz(),
                          uncorrJet.energy()
                            );
                 */


                /*
                   unCorrJetPx_.push_back(uncorrJet.px());
                   unCorrJetPy_.push_back(uncorrJet.py());
                   unCorrJetPz_.push_back(uncorrJet.pz());
                   unCorrJetE_.push_back(uncorrJet.energy());
                 */
                jetArea_.push_back(jet->jetArea());


                // if reading text files, set jet 4-momentum
                // make correction using jecText files

                if(useJECText_) {
                        jecText_->setJetEta( uncorrJet.eta() );
                        jecText_->setJetPt ( uncorrJet.pt() );
                        jecText_->setJetE  ( uncorrJet.energy() );
                        jecText_->setJetA  ( jet->jetArea() );
                        jecText_->setRho   ( *(h_rho.product()) );
                        jecText_->setNPV   ( h_pv->size() );
                        Float_t corr_jet = jecText_->getCorrection();

                        /*
                           new( (*jetP4_)[nJet_-1]) TLorentzVector(uncorrJet.px()*corr_jet,
                                  uncorrJet.py()*corr_jet,
                                  uncorrJet.pz()*corr_jet,
                                  uncorrJet.energy()*corr_jet);
                         */
                        jetPx_.push_back(uncorrJet.px()*corr_jet);
                        jetPy_.push_back(uncorrJet.py()*corr_jet);
                        jetPz_.push_back(uncorrJet.pz()*corr_jet);
                        jetE_.push_back(uncorrJet.energy()*corr_jet);

                        jecUncText_->setJetEta( uncorrJet.eta() );
                        jecUncText_->setJetPt( corr_jet * uncorrJet.pt() );
                        jetCorrUncUp_.push_back(jecUncText_->getUncertainty(true));

                        jecUncText_->setJetEta( uncorrJet.eta() );
                        jecUncText_->setJetPt( corr_jet * uncorrJet.pt() );
                        jetCorrUncDown_.push_back(jecUncText_->getUncertainty(false));

                        // add px, py, pz , energy
                        jetPx_.push_back(uncorrJet.px()*corr_jet);
                        jetPy_.push_back(uncorrJet.py()*corr_jet);
                        jetPz_.push_back(uncorrJet.pz()*corr_jet);
                        jetE_.push_back(uncorrJet.energy()*corr_jet);


                }
                else{
                        /*
                           new( (*jetP4_)[nJet_-1]) TLorentzVector(jet->p4().px(),
                                  jet->p4().py(),
                                  jet->p4().pz(),
                                  jet->p4().energy());
                         */
                        jetPx_.push_back(jet->p4().px());
                        jetPy_.push_back(jet->p4().py());
                        jetPz_.push_back(jet->p4().pz());
                        jetE_.push_back(jet->p4().energy());
                }
                // get jet energy scale uncertainty and related input variables
                // fat jet uncertainty does not exist yet, if using database
                if(!useJECText_) {

                        jecUnc_->setJetEta(jet->eta());
                        jecUnc_->setJetPt(jet->pt());
                        jetCorrUncUp_.push_back(jecUnc_->getUncertainty(true));

                        jecUnc_->setJetEta(jet->eta());
                        jecUnc_->setJetPt(jet->pt());
                        jetCorrUncDown_.push_back(jecUnc_->getUncertainty(false));

                        // source of Uncertainty

                        const int nsrc = 11;
                        std::vector<float> temp_uncer;
                        temp_uncer.clear();
                        for (int isrc = 0; isrc < nsrc; isrc++) {

                                unc = vsrc[isrc];
                                unc->setJetPt(jet->pt());
                                unc->setJetEta(jet->eta());
                                float value = unc->getUncertainty(true);
                                temp_uncer.push_back(value);
                        };

                        total->setJetPt(jet->pt());
                        total->setJetEta(jet->eta());
                        float uncert = total->getUncertainty(true);
                        total_.push_back(uncert);
                        uncerSources_.push_back(temp_uncer);
                }


                jetCharge_.push_back(jet->charge());
                jetPartonFlavor_.push_back(jet->partonFlavour());
                jetHadronFlavor_.push_back(jet->hadronFlavour());

                if(isTHINJet_) {
                        float jpumva=0.;
                        jpumva= jet->userFloat("pileupJetId:fullDiscriminant");
                        PUJetID_.push_back(jpumva);
                        //  b  jet regrssion correction
                        bRegNNCorr_.push_back(jet->userFloat("bRegNNCorr"));
                        bRegNNResolution_.push_back(jet->userFloat("bRegNNResolution"));

                        isPUJetIDLoose_.push_back(bool(jet->userInt("pileupJetId:fullId") & (1 << 2)));
                        isPUJetIDMedium_.push_back(bool(jet->userInt("pileupJetId:fullId") & (1 << 1)));
                        isPUJetIDTight_.push_back(bool(jet->userInt("pileupJetId:fullId") & (1 << 0)));
                }

                if(isTHINJet_) {
                        if (runOn2016_) {
                                jetPassIDLoose_.push_back(bool(jet->userInt("looseJetID_2016")));
                        }
                        if (runOn2017_) {
                                jetPassIDTight_.push_back(bool(jet->userInt("tightJetID_2017")));
                        }
                        if (runOn2018_) {
                                jetPassIDTight_.push_back(bool(jet->userInt("tightJetID_2018")));
                        }
                        jetCEmEF_.push_back(jet->userFloat("CEMF_"));
                        jetCHadEF_.push_back(jet->userFloat("CHF_"));
                        jetNEmEF_.push_back(jet->userFloat("NEMF_"));
                        jetNHadEF_.push_back(jet->userFloat("NHF_"));
                        jetMuoEF_.push_back(jet->userFloat("MUF_"));
                        jetCMulti_.push_back(jet->userFloat("CHM_"));
                        jetNMultiplicity_.push_back(jet->userFloat("NumNeutralParticles_"));

                }
                if(!isTHINJet_) {
                        if (runOn2016_) {
                                std::map<std::string, bool> Pass = jet2017ID_.LooseJetCut_2016(*jet);
                                bool passOrNot = PassAll(Pass);
                                jetPassIDLoose_.push_back(passOrNot);
                        }
                        if (runOn2017_) {
                                std::map<std::string, bool> PassT = jet2017ID_.TightJetCut_2017(*jet);
                                bool passOrNotT = PassAll(PassT);
                                jetPassIDTight_.push_back(passOrNotT);
                        }
                        if (runOn2018_) {
                                std::map<std::string, bool> PassT = jet2017ID_.TightJetCut_2018(*jet);
                                bool passOrNotT = PassAll(PassT);
                                jetPassIDTight_.push_back(passOrNotT);
                        }

                        jetCEmEF_.push_back(jet->chargedEmEnergyFraction());
                        jetCHadEF_.push_back(jet->chargedHadronEnergyFraction());
                        jetNEmEF_.push_back(jet->neutralEmEnergyFraction());
                        jetMuoEF_.push_back(jet->muonEnergyFraction());
                        jetNHadEF_.push_back(jet->neutralHadronEnergyFraction());
                        jetCMulti_.push_back(jet->chargedMultiplicity());
                        jetNMultiplicity_.push_back(jet->neutralMultiplicity());
                }

                jetPhoEF_.push_back(jet->photonEnergyFraction());
                jetEleEF_.push_back(jet->electronEnergyFraction());
                jetChMuEF_.push_back(jet->chargedMuEnergyFraction());
                jetHFHadEF_.push_back(jet->HFHadronEnergyFraction());
                jetHFEMEF_.push_back(jet->HFEMEnergyFraction());
                jetHOEnergy_.push_back(jet->hoEnergy());
                jetHOEF_.push_back(jet->hoEnergyFraction());
                jetEleMultiplicity_.push_back(jet->electronMultiplicity());
                jetMuoMultiplicity_.push_back(jet->muonMultiplicity());
                jetCHHadMultiplicity_.push_back(jet->chargedHadronMultiplicity());
                jetPhMultiplicity_.push_back(jet->photonMultiplicity());
                jetNHadMulplicity_.push_back(jet->neutralHadronMultiplicity());
                jetHFHadMultiplicity_.push_back(jet->HFHadronMultiplicity());
                jetHFEMMultiplicity_.push_back(jet->HFEMMultiplicity());


                // b-tagging

                jetSSV_.push_back(jet->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
                jetCSV_.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
                jetSSVHE_.push_back(jet->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
                jetCISVV2_.push_back(jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
                jetDeepCSV_b_.push_back(jet->bDiscriminator("pfDeepCSVJetTags:probb")+jet->bDiscriminator("pfDeepCSVJetTags:probbb"));
                jetDeepCSV_c_.push_back(jet->bDiscriminator("pfDeepCSVJetTags:probc"));
                jetDeepCSV_udsg_.push_back(jet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
                jetTCHP_.push_back(jet->bDiscriminator("pfTrackCountingHighPurBJetTags"));
                jetTCHE_.push_back(jet->bDiscriminator("pfTrackCountingHighEffBJetTags"));
                jetJP_.push_back(jet->bDiscriminator("pfJetProbabilityBJetTags"));
                jetJBP_.push_back(jet->bDiscriminator("pfJetBProbabilityBJetTags"));


                if(isCA15PuppiJet_) {

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
                        std::vector<edm::Ptr<reco::Candidate> > constituentPtrs = jet->getJetConstituents();

                        // if (!minimal && data->size()<2)   commented by raman becuase this was not needed here,
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

                                if (alljets.size()>0) {
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

                                }

                        } // if not minimal and fewer than 2
                          // End of ECFs computation


                        // double b-tagger calculation for the CA15 jets
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
                        float massPruned_                    = jet->p4().mass();//jet->m;
                        float flavour_                       = -1;//j.partonFlavor();   // they're spectator variables
                        float nbHadrons_                     = -1;//j.hadronFlavor(); //
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
                        ECF_2_3_10.push_back(p_jet->get_ecf(2,3,1));
                        ECF_1_2_10.push_back(p_jet->get_ecf(1,2,1));

                        jetTau1_.push_back(jet->userFloat("NjettinessCA15Puppi:tau1"));
                        jetTau2_.push_back(jet->userFloat("NjettinessCA15Puppi:tau2"));
                        jetTau3_.push_back(jet->userFloat("NjettinessCA15Puppi:tau3"));
                        jetTau4_.push_back(jet->userFloat("NjettinessCA15Puppi:tau2"));
                        jetSDmass_.push_back(jet->userFloat("ca15PFJetsPuppiSoftDropMass"));


                }//     if(isCA15PuppiJet_){


                // This is the main jet collection which will be used for the monoH analysis
                if(isFATJet_) {

                        // CHS subjettiness
                        jetCHSTau1_.push_back(jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));
                        jetCHSTau2_.push_back(jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2"));
                        jetCHSTau3_.push_back(jet->userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));


                        TLorentzVector temp_CHS;
                        temp_CHS.SetPtEtaPhiM(jet->userFloat("ak8PFJetsCHSValueMap:pt"),
                                              jet->userFloat("ak8PFJetsCHSValueMap:eta"),
                                              jet->userFloat("ak8PFJetsCHSValueMap:phi"),
                                              jet->userFloat("ak8PFJetsCHSValueMap:mass"));
                        //new( (*jetCHSP4_)[nJet_-1]) TLorentzVector(temp_CHS);


                        jetCHSPx_.push_back(temp_CHS.Px());
                        jetCHSPy_.push_back(temp_CHS.Py());
                        jetCHSPz_.push_back(temp_CHS.Pz());
                        jetCHSE_.push_back(temp_CHS.E());

                        //Puppi subjettiness
                        jetTau1_.push_back(jet->userFloat("NjettinessAK8Puppi:tau1"));
                        jetTau2_.push_back(jet->userFloat("NjettinessAK8Puppi:tau2"));
                        jetTau3_.push_back(jet->userFloat("NjettinessAK8Puppi:tau3"));
                        jetTau4_.push_back(jet->userFloat("NjettinessAK8Puppi:tau4"));
                        //jetSDmass_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropMass"));

                        // deep DoubleB tagger

                        jet_probQCDb_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probQCD"));
                        jet_probHbb_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb"));
                        jet_probQCDc_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probQCD"));
                        jet_probHcc_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleCvLJetTags:probHcc"));
                        jet_probHbbc_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHbb"));
                        jet_probHccb_.push_back(jet->bDiscriminator("pfMassIndependentDeepDoubleCvBJetTags:probHcc"));

                        jet_prob_bbvsLight_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:bbvsLight"));
                        jet_prob_ccvsLight_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ccvsLight"));
                        jet_prob_TvsQCD_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:TvsQCD"));
                        jet_prob_ZHccvsQCD_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHccvsQCD"));
                        jet_prob_WvsQCD_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:WvsQCD"));
                        jet_prob_ZHbbvsQCD_.push_back(jet->bDiscriminator("pfMassDecorrelatedDeepBoostedDiscriminatorsJetTags:ZHbbvsQCD"));


                        //jet__.push_back(jet->bDiscriminator(""));
//
                        N2_Beta1_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN2"));
                        N3_Beta1_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb1AK8PuppiSoftDropN3"));
                        N2_Beta2_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN2"));
                        N3_Beta2_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropValueMap:nb2AK8PuppiSoftDropN3"));




                        //   using a different way to get corrected pruned/softdrop mass
                        // if reading global tag
                        float corr=-1;

                        if(!useJECText_) {
                                // pruned mass: CHS
                                std::vector<pat::Jet>::const_iterator jetForPrunedMass =
                                        find_if(jetsForPrunedMass.begin(),
                                                jetsForPrunedMass.end(),
                                                [&jet](const pat::Jet& item)->bool {
                                        return fabs(jet->correctedP4(0).pt()-item.correctedP4(0).pt())<1e-3;
                                });

                                if(jetForPrunedMass!=jetsForPrunedMass.end())
                                        corr = jetForPrunedMass->pt()/jetForPrunedMass->correctedP4(0).pt();
                        }

                        else if(useJECText_) {
                                // pruned mass: CHS
                                prunedjecText_->setJetEta( uncorrJet.eta() );
                                prunedjecText_->setJetPt ( uncorrJet.pt() );
                                prunedjecText_->setJetE  ( uncorrJet.energy() );
                                prunedjecText_->setJetA  ( jet->jetArea() );
                                prunedjecText_->setRho   ( *(h_rho.product()) );
                                prunedjecText_->setNPV   ( h_pv->size() );
                                corr = prunedjecText_->getCorrection();

                        }

                        if(corr<0) {
                                jetCHSPRmassL2L3Corr_.push_back(DUMMY);
                                jetCHSSDmassL2L3Corr_.push_back(DUMMY);
                        }
                        else{
                                jetCHSPRmassL2L3Corr_.push_back(corr*jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass"));
                                jetCHSSDmassL2L3Corr_.push_back(corr*jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass"));
                        }

                        jetCHSSDmass_.push_back(jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass"));
                        jetCHSPRmass_.push_back(jet->userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass"));


                } // only for FAT jets in default MiniAOD (AK8Puppi)


                // if this is a AK8 jet
                if(!isTHINJet_ && !isAK4PuppiJet_) {
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

                                for(unsigned int n_2ndvtx=0; n_2ndvtx< nSV; n_2ndvtx++)
                                        jet_SVMass_float.push_back(candSVTagInfo->secondaryVertex(n_2ndvtx).p4().mass());

                        } // if there is tagging information
                        if(nSV==0)
                                jet_SVMass_float.push_back(DUMMY);
                        jet_nSV_.push_back(nSV);
                        jet_SVMass_.push_back(jet_SVMass_float);

                } /// if its FATjet or AK8PuppiJet or CA15PuppiJet

                // softdrop subjets
                if(isFATJet_ || isAK8PuppiJet_ || isCA15PuppiJet_) {
                        pat::Jet const *jetptr = &*jet;

                        //	std::cout<<" working before SoftDrop "<<std::endl;

                        auto wSubjets = isFATJet_? jetptr->subjets("SoftDropPuppi"): jetptr->subjets("SoftDrop");
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
                        TLorentzVector softdrop_raw_l4;
                        softdrop_raw_l4.SetPxPyPzE(0,0,0,0);

                        for ( auto const & iw : wSubjets )
                        {

                                nSubSoftDropjets++;

                                // build genjet softdrop mass
                                if (iw->genJet()) {
                                        genjet_softdrop_l4 +=TLorentzVector(
                                                iw->genJet()->p4().px(),
                                                iw->genJet()->p4().py(),
                                                iw->genJet()->p4().pz(),
                                                iw->genJet()->p4().energy()
                                                );
                                }

                                softdrop_raw_l4 += TLorentzVector(iw->correctedP4(0).px(),
                                                                  iw->correctedP4(0).py(),
                                                                  iw->correctedP4(0).pz(),
                                                                  iw->correctedP4(0).energy());

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

                        new( (*jetSDRawP4_)[nJet_-1]) TLorentzVector(softdrop_raw_l4);

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

                        if(isFATJet_)
                                //jetSDmass_.push_back(nSubSoftDropjets==0? DUMMY:softdrop_raw_l4.M());
                                jetSDmass_.push_back(jet->userFloat("ak8PFJetsPuppiSoftDropMass"));
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

bool jet_extra = false;

void
jetTree::SetBranches(){

        AddBranch(&nJet_,   "nJet");
        //AddBranch(&jetP4_,       "jetP4");
        AddBranch(&jetPx_, "jetPx");
        AddBranch(&jetPy_, "jetPy");
        AddBranch(&jetPz_, "jetPz");
        AddBranch(&jetE_, "jetEnergy");

        AddBranch(&jetRho_, "jetRho");
        AddBranch(&jetNPV_, "jetNPV");
        AddBranch(&uncerSources_,"jetUncSources");
        AddBranch(&total_,"jetUncTotal");

        AddBranch(&jetCEmEF_,  "jetCEmEF");
        AddBranch(&jetCHadEF_, "jetCHadEF");
        AddBranch(&jetNEmEF_,  "jetNEmEF");
        AddBranch(&jetNHadEF_, "jetNHadEF");
        AddBranch(&jetMuoEF_,  "jetMuoEF");
        AddBranch(&jetCMulti_, "jetCMulti");
        AddBranch(&jetNMultiplicity_,"jetNMultiplicity");
        AddBranch(&genjetpx_,"genjetpx");
        AddBranch(&genjetpy_,"genjetpy");
        AddBranch(&genjetpz_,"genjetpz");
        AddBranch(&genjetE_,"genjetE");

        if(jet_extra) {
                //AddBranch(&jetP4_,       "jetP4");
                //AddBranch(&genjetP4_,   "genjetP4"); // this is no longer needed as individual component is already there,
                AddBranch(&genjetEM_,  "genjetEM");
                AddBranch(&genjetHAD_, "genjetHAD");
                AddBranch(&genjetINV_, "genjetINV");
                AddBranch(&genjetAUX_, "genjetAUX");
                AddBranch(&matchedDR_, "matchedDR");
                AddBranch(&jetRawFactor_, "jetRawFactor");
                //AddBranch(&unCorrJetP4_, "unCorrJetP4");
                AddBranch(&unCorrJetPx_, "unCorrJetPx");
                AddBranch(&unCorrJetPy_, "unCorrJetPy");
                AddBranch(&unCorrJetPz_, "unCorrJetPz");
                AddBranch(&unCorrJetE_, "unCorrJetE");
                AddBranch(&jetArea_,        "jetArea");
                AddBranch(&jetCharge_,       "jetCharge");
                AddBranch(&jetPartonFlavor_, "jetPartonFlavor");
                AddBranch(&jetSSV_,   "jetSSV");
                AddBranch(&jetCSV_,   "jetCSV");
                AddBranch(&jetSSVHE_, "jetSSVHE");
                AddBranch(&jetTCHP_,  "jetTCHP");
                AddBranch(&jetTCHE_,  "jetTCHE");
                AddBranch(&jetJP_,    "jetJP");
                AddBranch(&jetJBP_,   "jetJBP");
                AddBranch(&jetEleMultiplicity_,"jetEleMulti");
                AddBranch(&jetMuoMultiplicity_,"jetMuoMulti");
                AddBranch(&jetCHHadMultiplicity_,"jetCHHadMultiplicity");
                AddBranch(&jetPhMultiplicity_,"jetPhMultiplicity");
                AddBranch(&jetNHadMulplicity_,"jetNHadMulplicity");
                AddBranch(&jetPhoEF_,  "jetPhoEF");
                AddBranch(&jetEleEF_,  "jetEleEF");
        }
        AddBranch(&jetCorrUncUp_,   "jetCorrUncUp");
        AddBranch(&jetCorrUncDown_, "jetCorrUncDown");
        AddBranch(&jetHadronFlavor_, "jetHadronFlavor");

        if (runOn2018_) {
                AddBranch(&jetPassIDTight_,  "jetPassIDTight");
        }
        else if (runOn2017_) {
                AddBranch(&jetPassIDTight_,  "jetPassIDTight");
        }
        else if (runOn2016_) {
                AddBranch(&jetPassIDLoose_,  "jetPassIDLoose");
        }

        if(isTHINJet_) {
                AddBranch(&jetCISVV2_,"jetCISVV2");
                AddBranch(&jetDeepCSV_b_,"jetDeepCSV_b");
                AddBranch(&jetDeepCSV_c_,"jetDeepCSV_c");
                AddBranch(&jetDeepCSV_udsg_,"jetDeepCSV_udsg");
                AddBranch(&PUJetID_,   "PUJetID");
                AddBranch(&isPUJetIDLoose_,  "isPUJetIDLoose");
                AddBranch(&isPUJetIDMedium_, "isPUJetIDMedium");
                AddBranch(&isPUJetIDTight_,  "isPUJetIDTight");
                AddBranch(&bRegNNCorr_,"bRegNNCorr");
                AddBranch(&bRegNNResolution_,"bRegNNResolution");
        }

        if(isFATJet_ || isAK8PuppiJet_ || isCA15PuppiJet_) {
                if (jet_extra) {

                        AddBranch(&jet_nSV_,     "jet_nSV");
                        AddBranch(&jet_SVMass_,  "jet_SVMass");

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
                }

                if (isFATJet_) {
                        AddBranch(&jet_DoubleSV_,"jet_DoubleSV");
                        AddBranch(&jet_probQCDb_,"jet_probQCDb");
                        AddBranch(&jet_probHbb_,"jet_probHbb");
                        AddBranch(&jet_probQCDc_,"jet_probQCDc");
                        AddBranch(&jet_probHcc_,"jet_probHcc");
                        AddBranch(&jet_probHbbc_,"jet_probHbbc");
                        AddBranch(&jet_probHccb_,"jet_probHccb");


                        AddBranch(&jet_prob_bbvsLight_, "jet_prob_bbvsLight");
                        AddBranch(&jet_prob_ccvsLight_, "jet_prob_ccvsLight");
                        AddBranch(&jet_prob_TvsQCD_, "jet_prob_TvsQCD");
                        AddBranch(&jet_prob_ZHccvsQCD_, "jet_prob_ZHccvsQCD");
                        AddBranch(&jet_prob_WvsQCD_, "jet_prob_WvsQCD");
                        AddBranch(&jet_prob_ZHbbvsQCD_, "jet_prob_ZHbbvsQCD");
                }

                AddBranch(&jetSDRawP4_, "jetSDRawP4");
                AddBranch(&jetSDmass_, "jetSDmass");

                AddBranch(&subjetSDHadronFlavor_, "subjetSDHadronFlavor");
                AddBranch(&subjetSDCSV_, "subjetSDCSV");

                if(isCA15PuppiJet_) {
                        AddBranch(&ca15_doublebtag, "_doublebtag");
                        AddBranch(&ECF_2_3_10, "ECF_2_3_10");
                        AddBranch(&ECF_1_2_10, "ECF_1_2_10");
                }

        } // for AK8 and CA15 jets

        if(isFATJet_) {
                AddBranch(&jetCHSSDmass_,         "jetCHSSDmass");
                AddBranch(&jetCHSPRmass_,         "jetCHSPRmass");
                AddBranch(&jetCHSPRmassL2L3Corr_, "jetCHSPRmassL2L3Corr");
                AddBranch(&jetCHSSDmassL2L3Corr_, "jetCHSSDmassL2L3Corr"); // newly added on 4.08.2019


                AddBranch(&jetCHSTau1_,   "jetCHSTau1");
                AddBranch(&jetCHSTau2_,   "jetCHSTau2");
                AddBranch(&jetCHSTau3_,   "jetCHSTau3");
                //AddBranch(&jetCHSP4_, "jetCHSP4");
                AddBranch(&jetCHSPx_, "jetCHSPx");
                AddBranch(&jetCHSPy_, "jetCHSPy");
                AddBranch(&jetCHSPz_, "jetCHSPz");
                AddBranch(&jetCHSE_, "jetCHSE");

                AddBranch(&jetTau1_,  "jetTau1");
                AddBranch(&jetTau2_,  "jetTau2");
                AddBranch(&jetTau3_,  "jetTau3");
                AddBranch(&jetTau4_, "jetTau4");

                AddBranch(&N2_Beta1_,"N2_Beta1_");
                AddBranch(&N3_Beta1_,"N3_Beta1_");
                AddBranch(&N2_Beta2_,"N2_Beta2_");
                AddBranch(&N3_Beta2_,"N3_Beta2_");


        } // only for AK8CHS jets


}


void
jetTree::Clear(){


        nJet_ =0;
        jetRho_=0;
        jetNPV_=0;

        //genjetP4_->Clear();
        genjetpx_.clear();
        genjetpy_.clear();
        genjetpz_.clear();
        genjetE_.clear();
        genjetEM_.clear();
        genjetHAD_.clear();
        genjetINV_.clear();
        genjetAUX_.clear();
        matchedDR_.clear();

        jetRawFactor_.clear();

        //jetP4_->Clear();
        //unCorrJetP4_->Clear();
        unCorrJetPx_.clear();
        unCorrJetPy_.clear();
        unCorrJetPz_.clear();
        unCorrJetE_.clear();

        jetPx_.clear();
        jetPy_.clear();
        jetPz_.clear();
        jetE_.clear();
        uncerSources_.clear();
        total_.clear();

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
        bRegNNCorr_.clear();
        bRegNNResolution_.clear();
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
        jetSSVHE_.clear();
        jetCISVV2_.clear();
        jetDeepCSV_b_.clear();
        jetDeepCSV_c_.clear();
        jetDeepCSV_udsg_.clear();
        jetTCHP_.clear();
        jetTCHE_.clear();
        jetJP_.clear();
        jetJBP_.clear();


        jetTau1_.clear();
        jetTau2_.clear();
        jetTau3_.clear();
        jetTau4_.clear();

        N2_Beta1_.clear();
        N3_Beta1_.clear();
        N2_Beta2_.clear();
        N3_Beta2_.clear();


        jetSDmass_.clear();
        jetSDRawP4_->Clear();

        // CHS related stuff

        jetCHSSDmass_.clear();
        jetCHSPRmass_.clear();
        jetCHSPRmassL2L3Corr_.clear();
        jetCHSSDmassL2L3Corr_.clear();


        jetCHSTau1_.clear();
        jetCHSTau2_.clear();
        jetCHSTau3_.clear();
        //jetCHSP4_->Clear();
        jetCHSPx_.clear();
        jetCHSPy_.clear();
        jetCHSPz_.clear();
        jetCHSE_.clear();

        // CA15 and ECFs
        ca15_doublebtag.clear();
        ECF_2_3_10.clear();
        ECF_1_2_10.clear();


        //jet  Hbb tagger for fat and add jet
        jet_probQCDb_.clear();
        jet_probHbb_.clear();
        jet_probQCDc_.clear();
        jet_probHcc_.clear();
        jet_probHbbc_.clear();
        jet_probHccb_.clear();

        jet_prob_bbvsLight_.clear();
        jet_prob_ccvsLight_.clear();
        jet_prob_TvsQCD_.clear();
        jet_prob_ZHccvsQCD_.clear();
        jet_prob_WvsQCD_.clear();
        jet_prob_ZHbbvsQCD_.clear();


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
