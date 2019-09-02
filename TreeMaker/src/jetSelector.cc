#include "ExoPieElement/TreeMaker/interface/jetSelector.h"


jetSelector::jetSelector()
{}

std::map<std::string, bool>  jetSelector::LooseJetCut_2016(const pat::Jet& jet){

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


std::map<std::string, bool>  jetSelector::TightJetCut_2016(const pat::Jet& jet){

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

std::map<std::string, bool>  jetSelector::TightJetCut_2017(const pat::Jet& jet){

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
          cuts["eta24"]       =  (jet.chargedHadronEnergyFraction() > 0.0) &&  (jet.chargedMultiplicity() > 0. );
        }
  }
  else if ((fabs(eta) > 2.7) && (fabs(eta) <= 3.0)) {
      cuts["etax"]        = (fabs(eta) > 2.7) && (fabs(eta) <= 3.0);
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() > 0.02 && jet.neutralEmEnergyFraction() < 0.99;
      cuts["nMulti"]        = jet.neutralMultiplicity() >= 2;
  }
  else {
      cuts["etax"]        = fabs(eta) > 3.0;
      cuts["nEmHF"]       = jet.neutralEmEnergyFraction() < 0.90;
      cuts["nHadEF"]      = jet.neutralHadronEnergyFraction() > 0.02;
      cuts["nMulti"]      = jet.neutralMultiplicity() > 10;
  }

  return cuts;
}
