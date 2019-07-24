#include "ExoPieElement/TreeMaker/interface/utils.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include <bitset>
//#include "TLorentzVector.h"
//#include <vector>
//#include <string>
//#include <map>

/*
class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.pt() > j.pt());
  }
};
*/
TLorentzVector Part2LorVec(reco::Candidate& cand){
  TLorentzVector* l = new TLorentzVector();
  l->SetPtEtaPhiM(cand.pt(),cand.eta(),cand.phi(),cand.mass());
  return (*l);
}

//When Selectors are fully developed,this must move to baseSelector class.
bool PassAll(std::map<std::string, bool> cutrecd){
  std::map<std::string, bool>::iterator iter= cutrecd.begin();
  bool decision =true ;
  for(;iter!=cutrecd.end();iter++){
    //    std::cout<<"-->"<<iter->first<<"\t"<<iter->second<<std::endl;	   
    decision = decision&&iter->second;     
  }
  return decision;
}


bool PassAllBut(std::string tag, std::map<std::string, bool> cutrecd){
  std::map<std::string, bool>::iterator iter= cutrecd.begin();
  bool decision =true ;
  for(;iter!=cutrecd.end();iter++){
    if(iter->first==tag)continue;
    decision = decision&&iter->second;
  }
 return decision;  
}



void getPFIsolation(double (&miniIso)[NISOPARS],
		    edm::Handle<pat::PackedCandidateCollection> pfcands,
		    reco::Candidate const* ptcl,  
		    const EffectiveAreas& eA_class, const double scEta,
		    const double rho,
		    const double r_iso_min, const double r_iso_max, const double kt_scale,
		    const bool charged_only) {

  if (ptcl->pt()<5.){
    for(int i=0; i<NISOPARS; i++)miniIso[i]=99999.;
    return;
  }
    double deadcone_nh(0.), deadcone_ch(0.), deadcone_ph(0.), deadcone_pu(0.), eA(0.);
    if(ptcl->isElectron()) {

      eA = eA_class.getEffectiveArea( scEta );
      if (fabs(ptcl->eta())>1.479) {deadcone_ch = 0.015; deadcone_pu = 0.015; deadcone_ph = 0.08;}

    } else if(ptcl->isMuon()) {

      eA = eA_class.getEffectiveArea( ptcl->eta() );
      deadcone_ch = 0.0001; deadcone_pu = 0.01; deadcone_ph = 0.01;deadcone_nh = 0.01;  

    } 


    double iso_nh(0.); double iso_ch(0.); 
    double iso_ph(0.); double iso_pu(0.);
    double ptThresh(0.5);
    if(ptcl->isElectron()) ptThresh = 0;
    double r_iso = TMath::Max(r_iso_min,TMath::Min(r_iso_max, kt_scale/ptcl->pt()));
    for (const pat::PackedCandidate &pfc : *pfcands) {
      if (abs(pfc.pdgId())<7) continue;

      double dr = deltaR(pfc, *ptcl);
      if (dr > r_iso) continue;
      
      //////////////////  NEUTRALS  /////////////////////////
      if (pfc.charge()==0){
        if (pfc.pt()>ptThresh) {
          /////////// PHOTONS ////////////
          if (abs(pfc.pdgId())==22) {
            if(dr < deadcone_ph) continue;
            iso_ph += pfc.pt();
	    /////////// NEUTRAL HADRONS ////////////
          } else if (abs(pfc.pdgId())==130) {
            if(dr < deadcone_nh) continue;
            iso_nh += pfc.pt();
          }
        }
        //////////////////  CHARGED from PV  /////////////////////////
      } else if (pfc.fromPV()>1){
        if (abs(pfc.pdgId())==211) {
          if(dr < deadcone_ch) continue;
          iso_ch += pfc.pt();
        }
        //////////////////  CHARGED from PU  /////////////////////////
      } else {
        if (pfc.pt()>ptThresh){
          if(dr < deadcone_pu) continue;
          iso_pu += pfc.pt();
        }
      }
    }

    double eA_miniIso = eA * r_iso*r_iso / 0.09;
    double iso_dBeta(0.);
    double iso_eArea(0.);
    if (charged_only){
      iso_dBeta = iso_eArea = iso_ch;
    } 
    else {
      iso_dBeta = iso_ch + TMath::Max(0., iso_nh + iso_ph - 0.5*iso_pu);
      iso_eArea = iso_ch + TMath::Max(0., iso_nh + iso_ph - rho*eA_miniIso);
    }
    miniIso[0] = iso_ch;
    miniIso[1] = iso_nh;
    miniIso[2] = iso_ph;
    miniIso[3] = iso_pu;
    miniIso[4] = r_iso;
    miniIso[5] = iso_dBeta/ptcl->pt();
    miniIso[6] = iso_eArea/ptcl->pt();
    return;
}








/// Auxiliary function to select muons in the context of high-pt muon
/// analysis. Could be extended to other analyses, of course.
bool CustisTrackerMuon (const reco::Muon* recoMu, const reco::Vertex& vertex) {
  
  //bool isGlobal = false;
  bool isTracker = false;
  //bool muonChamberHit = false;
  bool matchedStations = false;
  bool relativeError = false;
  bool dBCut = false;
  bool longiCut = false;
  bool pixelHit = false;
  bool trackerLayers = false;
  
  //isGlobal = recoMu->isGlobalMuon();
  isTracker = recoMu->isTrackerMuon();
  
  matchedStations = (recoMu->numberOfMatchedStations() > 1);
	
  const reco::TrackRef& bestTrackRef = recoMu->muonBestTrack();
  const reco::TrackRef& innerTrackRef = recoMu->innerTrack();
  const reco::TrackRef& trackRef = recoMu->track();
  
  relativeError = (bestTrackRef->ptError()/bestTrackRef->pt()) < 0.3;
  dBCut         = (fabs(bestTrackRef->dxy(vertex.position())) < 0.2);  
  longiCut      = (fabs(bestTrackRef->dz(vertex.position())) < 0.5);
  
  if(innerTrackRef.isNonnull())
    pixelHit = (innerTrackRef->hitPattern().numberOfValidPixelHits() > 0);
  
  if(trackRef.isNonnull())
    trackerLayers = (trackRef->hitPattern().trackerLayersWithMeasurement() > 5);
  
  bool passed = (isTracker and 
		 matchedStations and 
		 relativeError and 
		       dBCut and 
		 longiCut and 
		       pixelHit and 
		 trackerLayers);
  
  
  return passed;
}
