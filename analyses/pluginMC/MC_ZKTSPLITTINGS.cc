// -*- C++ -*-
#include "Rivet/Analyses/MC_KTSPLITTINGS_BASE.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZKTSPLITTINGS : public MC_KTSPLITTINGS_BASE {
  public:

    /// Default constructor
    MC_ZKTSPLITTINGS()
      : MC_KTSPLITTINGS_BASE("MC_ZKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      _dR = (getOption("SCHEME") == "BARE") ? 0.0 : 0.2;
      _lepton = (getOption("LMODE") == "MU") ? PID::MUON : PID::ELECTRON;

      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETALMAX", 3.5);
      const double ptcut = getOption<double>("PTLMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      DileptonFinder zfinder(91.2*GeV, _dR, cut && Cuts::abspid == _lepton, Cuts::massIn(65*GeV, 115*GeV));
      declare(zfinder, "DileptonFinder");

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);

      FastJets jetpro(zfinder.remainingFinalState(), JetAlg::KT, R);
      declare(jetpro, "Jets");

      MC_KTSPLITTINGS_BASE::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const DileptonFinder& zfinder = apply<DileptonFinder>(e, "DileptonFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;

      MC_KTSPLITTINGS_BASE::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_KTSPLITTINGS_BASE::finalize();
    }

    /// @}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    /// @{
    double _dR;
    PdgId _lepton;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_ZKTSPLITTINGS);

}
