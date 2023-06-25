// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "fastjet/contrib/SoftDrop.hh"

namespace Rivet {



  /// @brief MC validation analysis for jet events
  class MC_JETS : public MC_JetAnalysis {
  public:

    MC_JETS()
      : MC_JetAnalysis("MC_JETS", 4, "Jets")
    {    }


    void init() {
      // set ptcut from input option
      const double jetptcut = getOption<double>("PTJMIN", 20.0);
      _jetptcut = jetptcut * GeV;

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.4);

      // set clustering algorithm from input option
      FastJets::Algo clusterAlgo;
      const string algoopt = getOption("ALGO", "ANTIKT");

      if ( algoopt == "KT" ) {
	clusterAlgo = FastJets::KT;
      } else if ( algoopt == "CA" ) {
	clusterAlgo = FastJets::CA;
      } else if ( algoopt == "ANTIKT" ) {
	clusterAlgo = FastJets::ANTIKT;
      } else {
	MSG_WARNING("Unknown jet clustering algorithm option " + algoopt + ". "
		    "Defaulting to anti-kT");
	clusterAlgo = FastJets::ANTIKT;
      }
           
      FinalState fs;
      FastJets jetpro(fs, clusterAlgo, R);
      
      const string groomopt = getOption("GROOM", "");

      if (groomopt == "SD") {
        jetpro.addTrf(new fastjet::contrib::SoftDrop(0.0, 0.1));
      } else if (groomopt == "TRIM") {
        jetpro.addTrf(new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm, 0.2), fastjet::SelectorPtFractionMin(0.05)));
      } else if (groomopt != "") {
        MSG_WARNING("Unknown GROOM=" + groomopt + " option. Not applying jet grooming");
      }

      declare(jetpro, "Jets");
      MC_JetAnalysis::init();
    }


    void analyze(const Event& event) {
      MC_JetAnalysis::analyze(event);
    }


    void finalize() {
      MC_JetAnalysis::finalize();
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_JETS);

}
