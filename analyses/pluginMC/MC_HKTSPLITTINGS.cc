// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @brief MC validation analysis for higgs [-> tau tau] + jets events
  class MC_HKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_HKTSPLITTINGS()
      : MC_JetSplittings("MC_HKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETATAUMAX", 3.5);
      const double ptcut = getOption<double>("PTTAUMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      
      /// @todo Urk, abuse! Need explicit HiggsFinder and TauFinder
      ZFinder hfinder(FinalState(), cut, PID::TAU, 115*GeV, 135*GeV, 0.0, ZFinder::ClusterPhotons::NONE, ZFinder::AddPhotons::NO, 125*GeV);
      declare(hfinder, "Hfinder");

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(hfinder.remainingFinalState(), FastJets::KT, R);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = apply<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;
      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    //@}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_HKTSPLITTINGS);

}
