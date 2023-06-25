// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for W^+[enu]W^-[munu] + jets events
  class MC_WWKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_WWKTSPLITTINGS()
      : MC_JetSplittings("MC_WWKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      FinalState fs;

      // set FS cuts from input options
      const double etaecut = getOption<double>("ABSETAEMAX", 3.5);
      const double ptecut = getOption<double>("PTEMIN", 25.);

      Cut cute = Cuts::abseta < etaecut && Cuts::pT > ptecut*GeV;

      WFinder wenufinder(fs, cute, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wenufinder, "WenuFinder");

      VetoedFinalState wmnuinput;
      wmnuinput.addVetoOnThisFinalState(wenufinder);

      // set FS cuts from input options
      const double etamucut = getOption<double>("ABSETAMUMAX", 3.5);
      const double ptmucut = getOption<double>("PTMUMIN", 25.);

      Cut cutmu = Cuts::abseta < etamucut && Cuts::pT > ptmucut*GeV;
      
      WFinder wmnufinder(wmnuinput, cutmu, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wmnufinder, "WmnuFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(wenufinder)
          .addVetoOnThisFinalState(wmnufinder);

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);
      
      FastJets jetpro(jetinput, FastJets::KT, R);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const WFinder& wenufinder = apply<WFinder>(e, "WenuFinder");
      if (wenufinder.bosons().size()!=1) {
        vetoEvent;
      }

      const WFinder& wmnufinder = apply<WFinder>(e, "WmnuFinder");
      if (wmnufinder.bosons().size()!=1) {
        vetoEvent;
      }

      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    /// @}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_WWKTSPLITTINGS);

}
