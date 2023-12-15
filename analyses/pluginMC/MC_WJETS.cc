// -*- C++ -*-
#include "Rivet/Analyses/MC_JETS_BASE.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief MC validation analysis for W + jets events
  class MC_WJETS : public MC_JETS_BASE {
  public:

    /// Default constructor
    MC_WJETS(string name="MC_WJETS")
      : MC_JETS_BASE(name, 4, "Jets")
    {  }


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
      declare("MET", MissingMomentum());
      LeptonFinder lf(_dR, cut && Cuts::abspid == _lepton);
      declare(lf, "Leptons");

      // Set pT cut, jet alg, and clustering radius from input options
      _jetptcut = getOption<double>("PTJMIN", 20.0) * GeV;
      const double R = getOption<double>("R", 0.4);
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

      // Find jets with clustering radius from input option
      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(lf);
      FastJets fj(jetinput, clusterAlgo, R);
      declare(fj, "Jets");

      book(_h_W_jet1_deta, "W_jet1_deta", 50, -5.0, 5.0);
      book(_h_W_jet1_dR, "W_jet1_dR", 25, 0.5, 7.0);

      MC_JETS_BASE::init();
    }


    /// Do the analysis
    void analyze(const Event& event) {

      // MET cut
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      if (pmiss.pT() < 25*GeV) vetoEvent;

      // Identify the closest-matching l+MET to m == mW
      const Particles& ls = apply<LeptonFinder>(event, "Leptons").particles();
      const int ifound = closestMatchIndex(ls, pmiss, Kin::mass, 80.4*GeV, 60*GeV, 100*GeV);
      if (ifound < 0) vetoEvent;
      const Particle& l = ls[ifound];
      FourMomentum wmom = l.mom() + pmiss;

      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > _jetptcut);
      if (jets.size() > 0) {
        _h_W_jet1_deta->fill(wmom.eta()-jets[0].eta());
        _h_W_jet1_dR->fill(deltaR(wmom, jets[0].momentum()));
      }

      MC_JETS_BASE::analyze(event);
    }


    /// Finalize
    void finalize() {
      scale(_h_W_jet1_deta, crossSection()/picobarn/sumOfWeights());
      scale(_h_W_jet1_dR, crossSection()/picobarn/sumOfWeights());
      MC_JETS_BASE::finalize();
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
    Histo1DPtr _h_W_jet1_deta;
    Histo1DPtr _h_W_jet1_dR;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_WJETS);

}
