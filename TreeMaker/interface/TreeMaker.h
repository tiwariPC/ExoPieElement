#ifndef  TREE_MAKER_H
#define  TREE_MAKER_H

// user include files
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
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

#include "TTree.h"
#include "TFile.h"

class TreeMaker : public edm::EDAnalyzer {
   public:
      explicit TreeMaker(const edm::ParameterSet&);
      static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
      ~TreeMaker();


   private:
      virtual void beginJob() ;
      virtual void beginRun(const edm::Run&, const edm::EventSetup& ) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup& );
      virtual void endRun(const edm::Run&, const edm::EventSetup& );
      void beginLuminosityBlock(edm::LuminosityBlock const &, edm::EventSetup const &) override;
      edm::EDGetTokenT<GenLumiInfoHeader> genLumiHeaderToken_;
      void tokenize(std::string const &str, const char delim, std::vector<std::string> &out)
      {
         std::stringstream ss(str);
         std::string s;
         while (std::getline(ss, s, delim))
         {
            out.push_back(s);
         }
      }
      virtual void endJob() ;
      TFile* file;
      TTree* tree_;
      bool runOnSignal_;
      bool runOnrp2HDM_;
      bool runOnrpZpB_;
      std::string outFileName_ ;
      std::string scanId_;
      std::vector<std::string> mp_list;
      std::vector<std::string> chi_mass;
      std::vector<std::string> Zp_mass;
      int mass_A;
      int mass_a;
      int mass_chi;
      bool fillPUweightInfo_;
      bool fillEventInfo_;
      bool fillMetInfo_;
      bool fillTrigInfo_;
      bool fillFilterInfo_;

      bool fillGenInfo_;

      bool fillElecInfo_;
      bool fillMuonInfo_;
      bool fillTauInfo_;
      bool fillPhotInfo_;


      bool fillJetInfo_;
      bool fillFATJetInfo_;
      bool fillAK4PuppiJetInfo_;
      bool fillAK8PuppiJetInfo_;
      bool fillCA15PuppiJetInfo_;

      puweight *puweight_;
      eventInfo   *eventInfo_;
      patMetTree  *patMetTree_;

      patHltTree  *patHltTree_;
      patFilters  *patFilterTree_;

      genInfoTree *genInfoTree_;

      patElecTree *patElecTree_;
      patMuonTree *patMuTree_;
      hpstauInfo  *tauTree_;
      photonTree  *photonTree_;


      jetTree     *THINjetTree_;
      jetTree     *FATjetTree_;
      jetTree     *AK4PuppijetTree_;
      jetTree     *AK8PuppijetTree_;
      jetTree     *CA15PuppijetTree_;



};


#endif
