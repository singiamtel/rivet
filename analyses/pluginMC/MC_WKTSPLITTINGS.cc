// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @brief MC validation analysis for kt splitting scales in W + jets events
  class MC_WKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_WKTSPLITTINGS()
      : MC_JetSplittings("MC_WKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
		  _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETALMAX", 3.5);
      const double ptcut = getOption<double>("PTLMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      
      FinalState fs;
      WFinder wfinder(fs, cut, _lepton, 60.0*GeV, 100.0*GeV, 25.0*GeV, _dR);
      declare(wfinder, "WFinder");

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);
      
      FastJets jetpro(wfinder.remainingFinalState(), FastJets::KT, R);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }


    /// Do the analysis
    void analyze(const Event & e) {

      const WFinder& wfinder = apply<WFinder>(e, "WFinder");
      if (wfinder.bosons().size() != 1) {
        vetoEvent;
      }

      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    /// @}

  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    /// @{
    double _dR;
    PdgId _lepton;
    /// @}


  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_WKTSPLITTINGS);

}
