// -*- C++ -*-
#include "Rivet/Analyses/MC_KTSPLITTINGS_BASE.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// MC validation analysis for Higgs [-> tau tau] + jets events
  class MC_HKTSPLITTINGS : public MC_KTSPLITTINGS_BASE {
  public:

    /// Default constructor
    MC_HKTSPLITTINGS()
      : MC_KTSPLITTINGS_BASE("MC_HKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETATAUMAX", 3.5);
      const double ptcut = getOption<double>("PTTAUMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      /// @todo FS taus??
      DileptonFinder hfinder(125*GeV, 0.0, cut && Cuts::abspid == PID::TAU, Cuts::massIn(115*GeV, 135*GeV));

      declare(hfinder, "Hfinder");

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(hfinder.remainingFinalState(), JetAlg::KT, R);
      declare(jetpro, "Jets");

      MC_KTSPLITTINGS_BASE::init();
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const DileptonFinder& hfinder = apply<DileptonFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;
      MC_KTSPLITTINGS_BASE::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_KTSPLITTINGS_BASE::finalize();
    }

    /// @}

  };



  RIVET_DECLARE_PLUGIN(MC_HKTSPLITTINGS);

}
