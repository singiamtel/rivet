// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Higgs [-> tau tau] + jets events
  class MC_HJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_HJETS()
      : MC_JetAnalysis("MC_HJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETATAUMAX", 3.5);
      const double ptcut = getOption<double>("PTTAUMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;

      /// @todo Urk, abuse! Need explicit HiggsFinder (and TauFinder?)
      ZFinder hfinder(FinalState(), cut, PID::TAU, 115*GeV, 135*GeV,
		      0.0, PhotonOrigin::NONE, 125*GeV);
      declare(hfinder, "Hfinder");

      // set ptcut from input option
      const double jetptcut = getOption<double>("PTJMIN", 20.0);
      _jetptcut = jetptcut * GeV;

      // set clustering radius from input option
      const double R = getOption<double>("R", 0.4);

      // set clustering algorithm from input option
      JetAlg clusterAlgo;
      const string algoopt = getOption("ALGO", "ANTIKT");
      if ( algoopt == "KT" ) {
	clusterAlgo = JetAlg::KT;
      } else if ( algoopt == "CA" ) {
	clusterAlgo = JetAlg::CA;
      } else if ( algoopt == "ANTIKT" ) {
	clusterAlgo = JetAlg::ANTIKT;
      } else {
	MSG_WARNING("Unknown jet clustering algorithm option " + algoopt + ". Defaulting to anti-kT");
	clusterAlgo = JetAlg::ANTIKT;
      }

      FastJets jetpro(hfinder.remainingFinalState(), clusterAlgo, R);
      declare(jetpro, "Jets");

      book(_h_H_jet1_deta ,"H_jet1_deta", 50, -5.0, 5.0);
      book(_h_H_jet1_dR ,"H_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = apply<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;
      const double weight = 1.0;

      FourMomentum hmom(hfinder.bosons()[0].momentum());
      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size() > 0) {
        _h_H_jet1_deta->fill(hmom.eta()-jets[0].eta(), weight);
        _h_H_jet1_dR->fill(deltaR(hmom, jets[0].momentum()), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      normalize(_h_H_jet1_deta, crossSection()/picobarn);
      normalize(_h_H_jet1_dR, crossSection()/picobarn);
      MC_JetAnalysis::finalize();
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_H_jet1_deta;
    Histo1DPtr _h_H_jet1_dR;
    /// @}

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_HJETS);

}
