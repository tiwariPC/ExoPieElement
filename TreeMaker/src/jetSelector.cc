#include "ExoPieElement/TreeMaker/interface/jetSelector.h"


//jetSelector::jetSelector(const edm::ParameterSet& iConfig)
jetSelector::jetSelector()
{
//    runOn2017_(iConfig.getParameter<bool>("runOn2017")),
//    runOn2016_(iConfig.getParameter<bool>("runOn2016")),

}
runOn2017_ = True;
runOn2016_ = True;

std::map<std::string, bool>  jetSelector::MergedJetCut(const pat::Jet& jet){

  std::map<std::string, bool> cuts;
  float eta = jet.eta();

  //Kinematic and Fiducial Acceptance
  cuts["etax"]        = fabs(eta) < 2.5;
  cuts["muEF"]        = jet.muonEnergyFraction() < 0.99;
  cuts["phoEF"]       = jet.photonEnergyFraction() < 0.99;
  cuts["cEmEF"]       = jet.chargedEmEnergyFraction() < 0.99;
  cuts["nHadEF"]      = jet.neutralHadronEnergyFraction() < 0.99;
  cuts["cHadEF"]      = jet.chargedHadronEnergyFraction() > 0.0;


  // // not sure if this is needed
  cuts["trkSWBug"]    = !(fabs(eta)>1.0 && fabs(eta)<1.5 &&
    			  jet.chargedMultiplicity()/
			  (TMath::Max((float)0.1,(float)jet.neutralMultiplicity()))>2.0);
  return cuts;
}

std::map<std::string, bool>  jetSelector::LooseJetCut(const pat::Jet& jet){

  std::map<std::string, bool> cuts;
  float eta = jet.eta();

  //Kinematic and Fiducial Acceptance
  if (fabs(eta) < 2.7) {
      cuts["etax"]        = fabs(eta) <= 2.7;
      cuts["nHadEF"]      = jet.neutralHadronEnergyFraction() < 0.99;
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.99;
      cuts["cHadEF"]      = jet.chargedHadronEnergyFraction() > 0.0;
      cuts["cEmEF"]       = jet.chargedEmEnergyFraction() < 0.90;
      cuts["numConst"]    = (jet.chargedMultiplicity()+jet.neutralMultiplicity())>1;
      if (fabs(eta) < 2.4) {
          cuts["eta24"]       =  (jet.chargedHadronEnergyFraction() > 0.0) &&  (jet.chargedMultiplicity() > 0. ) &&  (jet.chargedEmEnergyFraction() < 0.99) ;
        }
  }
  else if ((fabs(eta) > 2.7) && (fabs(eta) <= 3.0)) {
      cuts["etax"]        = (fabs(eta) > 2.7) && (fabs(eta) <= 3.0);
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() > 0.01 && jet.neutralHadronEnergyFraction() < 0.98;
      cuts["nMulti"]        = jet.neutralMultiplicity() >= 2;
  }
  else {
      cuts["etax"]        = fabs(eta) > 3.0;
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.90;
      cuts["nMulti"]      = jet.neutralMultiplicity() > 10;
  }

  return cuts;
}



std::map<std::string, bool>  jetSelector::TightJetCut(const pat::Jet& jet){

  std::map<std::string, bool> cuts;
  float eta = jet.eta();

  //Kinematic and Fiducial Acceptance
  if (fabs(eta) < 2.7) {
      cuts["etax"]        = fabs(eta) <= 2.7;
      cuts["nHadEF"]      = jet.neutralHadronEnergyFraction() < 0.90;
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.90;
      cuts["cHadEF"]      = jet.chargedHadronEnergyFraction() > 0.0;
      cuts["cEmEF"]       = jet.chargedEmEnergyFraction() < 0.90;
      cuts["numConst"]    = (jet.chargedMultiplicity()+jet.neutralMultiplicity())>1;
      if (fabs(eta) < 2.4) {
          if (runOn2016_){
              cuts["eta24"]       =  (jet.chargedHadronEnergyFraction() > 0.0) &&  (jet.chargedMultiplicity() > 0. ) &&  (jet.chargedEmEnergyFraction() < 0.99) ;
          }
          if (runOn2017_){
              cuts["eta24"]       =  (jet.chargedHadronEnergyFraction() > 0.0) &&  (jet.chargedMultiplicity() > 0. );
          }

        }
  }
  else if ((fabs(eta) > 2.7) && (fabs(eta) <= 3.0)) {
      cuts["etax"]        = (fabs(eta) > 2.7) && (fabs(eta) <= 3.0);
      if (runOn2016_){
          cuts["nEmHF"]       = jet.neutralEmEnergyFraction() > 0.01 && jet.neutralHadronEnergyFraction() < 0.98;
      }
      if (runOn2017_){
          cuts["nEmHF"]       = jet.neutralEmEnergyFraction() > 0.02 && jet.neutralEmEnergyFraction() < 0.99;
      }
      cuts["nMulti"]        = jet.neutralMultiplicity() >= 2;
  }
  else {
      cuts["etax"]        = fabs(eta) > 3.0;
      if (runOn2016_){
          cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.90;
      }
      if (runOn2017_){
          cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.90;
          cuts["nHadEF"]      = jet.neutralHadronEnergyFraction() > 0.02;
      }
      cuts["nMulti"]      = jet.neutralMultiplicity() > 10;
  }

  return cuts;
}
