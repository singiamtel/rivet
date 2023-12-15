// -*- C++ -*-
#include "Rivet/Analyses/MC_KTSPLITTINGS_BASE.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// @brief MC validation analysis for e+mu W^+W^- + jets events
  class MC_WWKTSPLITTINGS : public MC_KTSPLITTINGS_BASE {
  public:

    /// Default constructor
    MC_WWKTSPLITTINGS()
      : MC_KTSPLITTINGS_BASE("MC_WWKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      declare("MET", MissingMomentum());

      // Find electrons with cuts from input options
      const double ETAECUT = getOption<double>("ABSETAEMAX", 3.5);
      const double PTECUT = getOption<double>("PTEMIN", 25.);
      const Cut cut_e = Cuts::abseta < ETAECUT && Cuts::pT > PTECUT*GeV;
      LeptonFinder ef(0.1, cut_e && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");

      // Find muons with cuts from input options
      const double ETAMUCUT = getOption<double>("ABSETAMUMAX", 3.5);
      const double PTMUCUT = getOption<double>("PTMUMIN", 25.);
      const Cut cut_m = Cuts::abseta < ETAMUCUT && Cuts::pT > PTMUCUT*GeV;
      LeptonFinder mf(0.2, cut_m && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      // Find jets with clustering radius from input option
      VetoedFinalState jetinput;
      jetinput
	.addVetoOnThisFinalState(ef)
	.addVetoOnThisFinalState(mf);
      const double R = getOption<double>("R", 0.6);
      FastJets fj(jetinput, JetAlg::KT, R);
      declare(fj, "Jets");

      MC_KTSPLITTINGS_BASE::init();
    }



    /// Do the analysis
    void analyze(const Event& event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      /// @note Dubious strategy, given there are two neutrinos...
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int iefound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      const Particles& mus = apply<LeptonFinder>(event, "Muons").particles();
      const int imfound = closestMatchIndex(mus, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);

      // Require two valid W candidates
      if (iefound < 0 || imfound < 0) vetoEvent;

      MC_KTSPLITTINGS_BASE::analyze(event);
    }


    /// Finalize
    void finalize() {
      MC_KTSPLITTINGS_BASE::finalize();
    }

    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_WWKTSPLITTINGS);

}
