// -*- C++ -*-
#include "Rivet/Analyses/MC_KTSPLITTINGS_BASE.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z[ee]Z[mumu] + jets events
  class MC_ZZKTSPLITTINGS : public MC_KTSPLITTINGS_BASE {
  public:

    /// Default constructor
    MC_ZZKTSPLITTINGS()
      : MC_KTSPLITTINGS_BASE("MC_ZZKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {

      // set FS cuts from input options
      const double etaecut = getOption<double>("ABSETAEMAX", 3.5);
      const double ptecut = getOption<double>("PTEMIN", 25.);
      Cut cute = Cuts::abseta < etaecut && Cuts::pT > ptecut*GeV;
      DileptonFinder zeefinder(91.2*GeV, 0.2, cute &&
                               Cuts::abspid == PID::ELECTRON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zeefinder, "ZeeFinder");

      VetoedFinalState zmminput;
      zmminput.vetoFinalState(zeefinder);

      // set FS cuts from input options
      const double etamucut = getOption<double>("ABSETAMUMAX", 3.5);
      const double ptmucut = getOption<double>("PTMUMIN", 25.);
      Cut cutmu = Cuts::abseta < etamucut && Cuts::pT > ptmucut*GeV;
      DileptonFinder zmmfinder(PromptFinalState(zmminput), 91.2*GeV, 0.2, cutmu &&
                               Cuts::abspid == PID::MUON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(zeefinder)
          .addVetoOnThisFinalState(zmmfinder);

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(jetinput, JetAlg::KT, R);
      declare(jetpro, "Jets");

      MC_KTSPLITTINGS_BASE::init();
    }


    /// Do the analysis
    void analyze(const Event & e) {
      const DileptonFinder& zeefinder = apply<DileptonFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size() != 1) vetoEvent;
      const DileptonFinder& zmmfinder = apply<DileptonFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size() != 1) vetoEvent;
      MC_KTSPLITTINGS_BASE::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_KTSPLITTINGS_BASE::finalize();
    }

    /// @}

  };



  RIVET_DECLARE_PLUGIN(MC_ZZKTSPLITTINGS);

}
