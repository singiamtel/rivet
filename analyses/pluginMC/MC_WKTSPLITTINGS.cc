// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
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

      // Use analysis options
      _dR = (getOption("SCHEME") == "BARE") ? 0.0 : 0.2;
      _lepton = (getOption("LMODE") == "MU") ? PID::MUON : PID::ELECTRON;
      const double ETACUT = getOption<double>("ABSETALMAX", 3.5);
      const double PTCUT = getOption<double>("PTLMIN", 25.);
      const Cut cut = Cuts::abseta < ETACUT && Cuts::pT > PTCUT*GeV;

      // Define projections
      LeptonFinder lf(_dR, cut && Cuts::abspid == _lepton);
      declare(lf, "Leptons");

      VetoedFinalState jetinput;
      jetinput.vetoFinalState(lf);
      const double R = getOption<double>("R", 0.6);
      FastJets fj(jetinput, JetAlg::KT, R);
      declare(fj, "Jets");

      MC_JetSplittings::init();
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& ls = apply<LeptonFinder>(event, "Leptons").particles();
      const int ifound = closestMassIndex(ls, pmiss, 80.4*GeV, 60*GeV, 100*GeV);
      if (ifound < 0) vetoEvent;

      MC_JetSplittings::analyze(event);
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


  RIVET_DECLARE_PLUGIN(MC_WKTSPLITTINGS);

}
