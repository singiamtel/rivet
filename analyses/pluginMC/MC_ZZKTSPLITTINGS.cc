// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z[ee]Z[mumu] + jets events
  class MC_ZZKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_ZZKTSPLITTINGS()
      : MC_JetSplittings("MC_ZZKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {

      // set FS cuts from input options
      const double etaecut = getOption<double>("ABSETAEMAX", 3.5);
      const double ptecut = getOption<double>("PTEMIN", 25.);

      Cut cute = Cuts::abseta < etaecut && Cuts::pT > ptecut*GeV;

      ZFinder zeefinder(FinalState(), cute, PID::ELECTRON, 65*GeV, 115*GeV,
                        0.2, PhotonOrigin::NODECAY, PhotonsAsConstituents::YES);
      declare(zeefinder, "ZeeFinder");

      VetoedFinalState zmminput;
      zmminput.addVetoOnThisFinalState(zeefinder);

      // set FS cuts from input options
      const double etamucut = getOption<double>("ABSETAMUMAX", 3.5);
      const double ptmucut = getOption<double>("PTMUMIN", 25.);

      Cut cutmu = Cuts::abseta < etamucut && Cuts::pT > ptmucut*GeV;
      
      ZFinder zmmfinder(zmminput, cutmu, PID::MUON, 65*GeV, 115*GeV,
                        0.2, PhotonOrigin::NODECAY, PhotonsAsConstituents::YES);
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(zeefinder)
          .addVetoOnThisFinalState(zmmfinder);

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);
      
      FastJets jetpro(jetinput, JetAlg::KT, R);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zeefinder = apply<ZFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size() != 1) vetoEvent;
      const ZFinder& zmmfinder = apply<ZFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size() != 1) vetoEvent;
      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    /// @}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_ZZKTSPLITTINGS);

}
