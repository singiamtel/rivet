// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief Demonstrate use of the Cutflow(s) classes
  class EXAMPLE_CUTFLOW : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(EXAMPLE_CUTFLOW);


    /// @name Analysis methods
    /// @{

    /// Set up projections and book histograms
    void init() {

      // Projections
      const DirectFinalState leps(Cuts::abseta < 4 && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      declare(leps, "Leptons");

      const FinalState fs(Cuts::abseta < 4);
      const FastJets jets(fs, JetAlg::Algo::ANTIKT, 0.4);
      declare(jets, "Jets");

      // Histograms
      // book(_s_cutflow, "cutflow");

      // Cut-flows
      book(_cutflows, {"JetLep", "DiBjet"}, {{"Jets", "Nlep", "pTlep1", "yLep1"},
                                             {"Jets", "0Lep", "2Jet", "Nbjets"}});
    }


    /// Do the analysis
    void analyze(const Event& event) {

      _cutflows->groupfillinit();

      const Jets jets = apply<JetFinder>(event, "Jets").jetsByPt();
      const Particles leps = apply<FinalState>(event, "Leptons").particlesByPt();

      if (jets.empty()) vetoEvent;
      _cutflows->groupfillnext();

      if (!leps.empty()) {
        _cutflows->fillnext("JetLep");
        if (!_cutflows->fillnext("JetLep", leps[0].pT() < 20*GeV)) vetoEvent;
        _cutflows->fillnext("JetLep", leps[0].absrap() < 2.5);
      } else {
        auto& mycf = _cutflows->binAt("DiBjet"s);
        mycf->fillnext();
        const Jets bjets = select(jets, hasBTag(Cuts::pT > 5*GeV));
        mycf->fillnext({jets.size() >= 2, bjets.size() == 2});
      }

    }

    /// Finalize
    void finalize() {
      MSG_INFO("Cut-flow:\n" << _cutflows);
    }

    /// @}


    /// Cut-flow counters
    CutflowsPtr _cutflows;

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(EXAMPLE_CUTFLOW);

}
