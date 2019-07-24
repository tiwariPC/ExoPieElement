// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Deepak Kumar
//         Created:  Wed, 10 Apr 2019 11:55:44 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "ExoPieElement/TreeMaker/interface/puweight.h"
#include "ExoPieElement/TreeMaker/interface/eventInfo.h"
#include "ExoPieElement/TreeMaker/interface/patMetTree.h"
#include "ExoPieElement/TreeMaker/interface/patHltTree.h"
#include "ExoPieElement/TreeMaker/interface/patFilters.h"
#include "ExoPieElement/TreeMaker/interface/genInfoTree.h"
#include "ExoPieElement/TreeMaker/interface/patElecTree.h"
#include "ExoPieElement/TreeMaker/interface/patMuonTree.h"
#include "ExoPieElement/TreeMaker/interface/hpstauInfo.h"
#include "ExoPieElement/TreeMaker/interface/photonTree.h"
#include "ExoPieElement/TreeMaker/interface/jetTree.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "ExoPieElement/TreeMaker/interface/jetTree.h"



#include "FWCore/ParameterSet/interface/ParameterSet.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class DemoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

    //  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
      


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      edm::EDGetTokenT<View<pat::Jet> > jetToken_;

      // ----------member data ---------------------------
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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
jetToken_( consumes<View<pat::Jet> >( iConfig.getParameter<InputTag> ( "JetTag" ) ) )
//jetToken_( iConfig.getParameter<InputTag> ( "JetTag" ) )
{
//produces<vector<pat::Jet> >(); 
  //now do what ever initialization is needed
   //usesResource("TFileService");
//jetToken_( consumes<View<pat::Jet> >( iConfig.getParameter<InputTag> ( "JetTag" ) ) )

}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

 Handle<View<pat::Jet> > jets;
 iEvent.getByToken( jetToken_, jets );

for(auto j = jets->begin(); j != jets->end(); ++j){
        if(j->pt() < 150.) continue;
        std::cout << "jetpT : " << j->pt() << std::endl;
        std::cout << "deepDoubleB = " << j->bDiscriminator("pfMassIndependentDeepDoubleBvLJetTags:probHbb") << std::endl;
}

//#ifdef THIS_IS_AN_EVENT_EXAMPLE
   //Handle<ExampleData> pIn;
   //iEvent.getByLabel("example",pIn);
//#endif
   
//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   //ESHandle<SetupData> pSetup;
   //iSetup.get<SetupRecord>().get(pSetup);
//#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
