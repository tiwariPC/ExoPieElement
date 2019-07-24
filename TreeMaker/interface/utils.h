#ifndef __UTILS_H_
#define __UTILS_H_
#include "TLorentzVector.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h" 
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include <vector>
#include <string>
#include <map>
class PtGreater {
  public:
  template <typename T> bool operator () (const T& i, const T& j) {
    return (i.pt() > j.pt());
  }
};


TLorentzVector Part2LorVec(reco::Candidate& cand);
bool PassAll(std::map<std::string, bool> cutrecd);
bool PassAllBut(std::string, std::map<std::string, bool>);

const int NISOPARS=7;
void getPFIsolation(double (&miniIso)[NISOPARS],
		    edm::Handle<pat::PackedCandidateCollection>,
		    reco::Candidate const*, 
		    const EffectiveAreas& eA_class,const double,
		    const double, const double, const double, const double, const bool);

bool CustisTrackerMuon (const reco::Muon* recoMu, const reco::Vertex& vertex) ;
#endif

