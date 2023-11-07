// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @brief MC validation analysis for jet events
  class MC_KTSPLITTINGS : public MC_JetSplittings {
  public:

    MC_KTSPLITTINGS()
      : MC_JetSplittings("MC_KTSPLITTINGS", 4, "Jets")
    {    }


  public:

    void init() {
      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(FinalState(), JetAlg::KT, R);
      declare(jetpro, "Jets");
      MC_JetSplittings::init();
    }


    void analyze(const Event& event) {
      MC_JetSplittings::analyze(event);
    }


    void finalize() {
      MC_JetSplittings::finalize();
    }

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_KTSPLITTINGS);

}
