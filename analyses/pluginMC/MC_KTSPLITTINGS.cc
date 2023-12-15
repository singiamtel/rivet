// -*- C++ -*-
#include "Rivet/Analyses/MC_KTSPLITTINGS_BASE.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @brief MC validation analysis for jet events
  class MC_KTSPLITTINGS : public MC_KTSPLITTINGS_BASE {
  public:

    MC_KTSPLITTINGS()
      : MC_KTSPLITTINGS_BASE("MC_KTSPLITTINGS", 4, "Jets")
    {    }


  public:

    void init() {
      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(FinalState(), JetAlg::KT, R);
      declare(jetpro, "Jets");
      MC_KTSPLITTINGS_BASE::init();
    }


    void analyze(const Event& event) {
      MC_KTSPLITTINGS_BASE::analyze(event);
    }


    void finalize() {
      MC_KTSPLITTINGS_BASE::finalize();
    }

  };


  RIVET_DECLARE_PLUGIN(MC_KTSPLITTINGS);

}
