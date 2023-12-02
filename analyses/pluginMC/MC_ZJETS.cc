// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z + jets events
  class MC_ZJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_ZJETS(string name = "MC_ZJETS")
      : MC_JetAnalysis(name, 4, "Jets")
	  {	  }


    /// @name Analysis methods
    /// @{

    /// Initialize
    void init() {
		  _dR=0.2;
      if (getOption("SCHEME") == "BARE")  _dR = 0.0;
		  _lepton=PID::ELECTRON;
      if (getOption("LMODE") == "MU")  _lepton = PID::MUON;

      // set FS cuts from input options
      const double etacut = getOption<double>("ABSETALMAX", 3.5);
      const double ptcut = getOption<double>("PTLMIN", 25.);

      Cut cut = Cuts::abseta < etacut && Cuts::pT > ptcut*GeV;
      DileptonFinder zfinder(91.2*GeV, _dR, cut && Cuts::abspid == _lepton, Cuts::massIn(66.0*GeV, 116.0*GeV));
      declare(zfinder, "DileptonFinder");

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

      FastJets jetpro(zfinder.remainingFinalState(), clusterAlgo, R);
      declare(jetpro, "Jets");

      book(_h_Z_jet1_deta ,"Z_jet1_deta", 50, -5, 5);
      book(_h_Z_jet1_dR ,"Z_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      MSG_TRACE("MC_ZJETS: running DileptonFinder");
      const DileptonFinder& zfinder = apply<DileptonFinder>(e, "DileptonFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const FourMomentum& zmom = zfinder.bosons()[0].momentum();
      MSG_TRACE("MC_ZJETS: have exactly one Z boson candidate");

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(Cuts::pT > _jetptcut);
      if (jets.size() > 0) {
        MSG_TRACE("MC_ZJETS: have at least one valid jet");
        _h_Z_jet1_deta->fill(zmom.eta()-jets[0].eta());
        _h_Z_jet1_dR->fill(deltaR(zmom, jets[0].momentum()));
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_Z_jet1_deta, crossSection()/picobarn/sumOfWeights());
      scale(_h_Z_jet1_dR, crossSection()/picobarn/sumOfWeights());
      MC_JetAnalysis::finalize();
    }

    /// @}


  protected:

    /// @name Parameters for specialised e/mu and dressed/bare subclassing
    /// @{
    double _dR;
    PdgId _lepton;
    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_Z_jet1_deta;
    Histo1DPtr _h_Z_jet1_dR;
    /// @}

  };

  // The hooks for the plugin system
  RIVET_DECLARE_PLUGIN(MC_ZJETS);
}
