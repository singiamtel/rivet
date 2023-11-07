// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_ZKTSPLITTINGS()
      : MC_JetSplittings("MC_ZKTSPLITTINGS", 4, "Jets")
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
      
      FinalState fs;
      ZFinder zfinder(fs, cut, _lepton, 65*GeV, 115*GeV, _dR,
		      PhotonOrigin::NODECAY, PhotonsAsConstituents::YES);
      declare(zfinder, "ZFinder");

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.6);
      
      FastJets jetpro(zfinder.remainingFinalState(), JetAlg::KT, R);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;

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


  RIVET_DECLARE_PLUGIN(MC_ZKTSPLITTINGS);
  
}
