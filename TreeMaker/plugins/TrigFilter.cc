// -*- C++ -*-
//
// Package:    TrigFilter
// Class:      TrigFilter
//
/**\class TrigFilter TrigFilter.cc ExoPieElement/TrigFilter/src/TrigFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Anil Pratap Singh,,,
//         Created:  Mon Sep 19 01:43:09 CEST 2011
// $Id: TrigFilter.cc,v 1.1 2012/03/15 09:58:20 lovedeep Exp $
//
//


// system include files
#include <memory>
#include<string>


// user include files

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/Common/interface/View.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//
// class declaration
namespace{
  template<typename T>
  edm::Handle<T> getHandle(const edm::Event& iEvent,const edm::EDGetTokenT<T>& token)
  {
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);
    return handle;
  }
}
//the functions which actually match the trigger objects and see if it passes
namespace{
  std::vector<const pat::TriggerObjectStandAlone*> getMatchedObjs(const float eta,const float phi,const std::vector<pat::TriggerObjectStandAlone>& trigObjs,const float maxDeltaR=0.1)
  {
    std::vector<const pat::TriggerObjectStandAlone*> matchedObjs;
    const float maxDR2 = maxDeltaR*maxDeltaR;
    for(auto& trigObj : trigObjs){
      const float dR2 = reco::deltaR2(eta,phi,trigObj.eta(),trigObj.phi());
      if(dR2<maxDR2) matchedObjs.push_back(&trigObj);
    }
    return matchedObjs;
  }
}

class TrigFilter : public edm::EDFilter {
   public:
      explicit TrigFilter(const edm::ParameterSet&);
      ~TrigFilter();
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      edm::EDGetTokenT<edm::TriggerResults>       trigResultToken;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection>     triggerObjectsToken;
      edm::EDGetTokenT<edm::View<pat::Electron> > elesToken_;
   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------
      edm::InputTag trigTag_;
      std::vector<std::string> triglist;
      bool isMC_;
      edm::InputTag patTrigObj_;

      bool isValidHltConfig_;
      HLTConfigProvider  hltConfigProvider_;
  TH1F* event_counter_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TrigFilter::TrigFilter(const edm::ParameterSet& iConfig):
 elesToken_(consumes<edm::View<pat::Electron> >(edm::InputTag("slimmedElectrons"))),
 trigTag_(iConfig.getParameter<edm::InputTag> ("TrigTag")),
 triglist(iConfig.getParameter<std::vector< std::string > >("TrigPaths")),
 isMC_( iConfig.getParameter<bool>( "isMC_" ) ),
 patTrigObj_(edm::InputTag("slimmedPatTrigger"))
{
   //now do what ever initialization is needed
   isValidHltConfig_ = false;
   trigResultToken = consumes<edm::TriggerResults>(trigTag_);
   triggerObjectsToken = consumes<std::vector<pat::TriggerObjectStandAlone> >(patTrigObj_);
   edm::Service<TFileService> fs;
   event_counter_ = fs->make<TH1F>("event_counter_", "event_counter_", 2,0.5,2.5);

}


TrigFilter::~TrigFilter()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
TrigFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  event_counter_->Fill(1);
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  iEvent.getByToken(triggerObjectsToken,triggerObjects);

  edm::Handle<edm::TriggerResults> trigResults;

  if (not iEvent.getByToken(trigResultToken, trigResults)) {
    std::cout << ">>> TRIGGER collection does not exist !!!\n";
    return false;
  }
  const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);
  bool decision = false;
  for (unsigned int i=0; i<trigResults->size(); i++){
    std::string trigName   = trigNames.triggerName(i);
    bool  trigResult = trigResults->accept(i);
    //std::vector<std::string>::const_iterator iter = triglist.begin();
    //Not interested in presclaed paths.
    bool prescaled = 0;
    const unsigned int prescaleSize=  hltConfigProvider_.prescaleSize();
    for(unsigned int ps=0; ps<prescaleSize;ps++){
      const unsigned int prescaleValue=hltConfigProvider_.prescaleValue(ps,trigName);
      if(prescaleValue !=1 )prescaled = true;
    }
    if(prescaled) continue;
    bool  ele_L1DoubleEG_trig = false;
    if (!isMC_){
      std::string trigName_1234 = trigName.substr(0, trigName.find("_v"));
      if( find(triglist.begin(), triglist.end(), trigName_1234)  == triglist.end() ) continue;
      std::string trig_ele_comp =  "HLT_Ele32_WPTight_Gsf_L1DoubleEG";
      //std::cout<<"trigName_1234:  "<<trigName_1234<<std::endl;
      if((trigName_1234.compare(trig_ele_comp)) == 0)
        {
            auto trigResultsHandle = getHandle(iEvent, trigResultToken);
            auto trigObjsHandle = getHandle(iEvent, triggerObjectsToken);
            auto elesHandle = getHandle(iEvent, elesToken_);
            // so the filter names are all packed in miniAOD so we need to create a new collection of them which are unpacked
            std::vector<pat::TriggerObjectStandAlone> unpackedTrigObjs;
            for (auto& trigObj : *trigObjsHandle)
            {
                unpackedTrigObjs.push_back(trigObj);
                unpackedTrigObjs.back().unpackFilterLabels(iEvent, *trigResultsHandle);
            }
            //std::cout << "checking eles " << std::endl;
            for (auto& ele : *elesHandle)
            {
                // the eta/phi of e/gamma trigger objects is the supercluster eta/phi
                const float eta = ele.superCluster()->eta();
                const float phi = ele.superCluster()->phi();

                // now match ALL objects in a cone of DR<0.1
                // it is important to match all objects as there are different ways to reconstruct the same electron
                // eg, L1 seeded, unseeded, as a jet etc and so you want to be sure you get all possible objects

                std::vector<const pat::TriggerObjectStandAlone*> matchedTrigObjs = getMatchedObjs(eta, phi, unpackedTrigObjs, 0.1);
                for (const auto trigObj : matchedTrigObjs)
                {
                    // now just check if it passes the two filters
                    if(trigObj->hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter") && trigObj->hasFilterLabel("hltEGL1SingleEGOrFilter") )
                    {
                        //std::cout << " ele " << ele.et() << " " << eta << " " << phi << " passes HLT_Ele32_WPTight_Gsf" << std::endl;
                        ele_L1DoubleEG_trig = true;
                    }
                }
            }
        }

      }
    if(ele_L1DoubleEG_trig) trigResult = trigResult && ele_L1DoubleEG_trig;
    decision = decision||trigResult;
  }
  return decision;
}

// ------------ method called once each job just before starting event loop  ------------
void
TrigFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
TrigFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool
TrigFilter::beginRun(edm::Run& r, edm::EventSetup const& iSetup)
{
  bool isConfigChanged = false;
  isValidHltConfig_ = hltConfigProvider_.init( r, iSetup, trigTag_.process(), isConfigChanged );
//  return isValidHltConfig_;
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool
TrigFilter::endRun(edm::Run&, edm::EventSetup const& iSetup)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool
TrigFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool
TrigFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TrigFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(TrigFilter);
