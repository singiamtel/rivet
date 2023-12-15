// -*- C++ -*-
#include "Rivet/Analyses/MC_JETS_BASE.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/LeptonFinder.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// @brief MC validation analysis for e+mu W^+W^- + jets events
  class MC_WWJETS : public MC_JETS_BASE {
  public:

    /// Default constructor
    MC_WWJETS()
      : MC_JETS_BASE("MC_WWJETS", 4, "Jets")
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
      LeptonFinder ef(0.2, cut_e && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");

      // Find muons with cuts from input options
      const double ETAMUCUT = getOption<double>("ABSETAMUMAX", 3.5);
      const double PTMUCUT = getOption<double>("PTMUMIN", 25.);
      const Cut cut_m = Cuts::abseta < ETAMUCUT && Cuts::pT > PTMUCUT*GeV;
      LeptonFinder mf(0.2, cut_m && Cuts::abspid == PID::MUON);
      declare(mf, "Muons");

      // Set pT cut, jet alg, and clustering radius from input options
      const double jetptcut = getOption<double>("PTJMIN", 20.0);
      _jetptcut = jetptcut * GeV;
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
      jetinput
      .addVetoOnThisFinalState(ef)
      .addVetoOnThisFinalState(mf);
      FastJets fj(jetinput, clusterAlgo, R);
      declare(fj, "Jets");


      // Histogram correlations with jets and global quantities
      book(_h_WW_jet1_deta ,"WW_jet1_deta", 70, -7.0, 7.0);
      book(_h_WW_jet1_dR ,"WW_jet1_dR", 25, 1.5, 7.0);
      book(_h_We_jet1_dR ,"We_jet1_dR", 25, 0.0, 7.0);
      //
      book(_h_HT ,"HT", logspace(100, 100.0, 0.5*(sqrtS()>0.?sqrtS():14000.)));
      book(_h_jets_m_12 ,"jets_m_12", logspace(100, 1.0, 0.25*(sqrtS()>0.?sqrtS():14000.)));

      MC_JETS_BASE::init();
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

      // Get momenta
      FourMomentum pe = es[iefound].mom();
      FourMomentum pm = mus[imfound].mom();
      FourMomentum pww = pe + pm + pmiss;

      // Get jets and compute obs relative to the WW objects
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > _jetptcut);
      if (jets.size() > 0) {
        _h_WW_jet1_deta->fill(pww.eta()-jets[0].eta());
        _h_WW_jet1_dR->fill(deltaR(pww, jets[0].momentum()));
        _h_We_jet1_dR->fill(deltaR(pe, jets[0].momentum()));
      }
      const double HT = sum(jets, Kin::pT, pe.pT() + pm.pT() + pmiss.pT());
      if (HT > 0.0) _h_HT->fill(HT/GeV);

      if (jets.size() > 1) {
        const FourMomentum jet1 = jets[0].momentum();
        const FourMomentum jet2 = jets[1].momentum();
        _h_jets_m_12->fill((jet1+jet2).mass()/GeV);
      }

      MC_JETS_BASE::analyze(event);
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/picobarn/sumOfWeights();
      scale(_h_WW_jet1_deta, norm);
      scale(_h_WW_jet1_dR, norm);
      scale(_h_We_jet1_dR, norm);
      scale(_h_jets_m_12, norm);
      scale(_h_HT, norm);
      MC_JETS_BASE::finalize();
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_WW_jet1_deta;
    Histo1DPtr _h_WW_jet1_dR;
    Histo1DPtr _h_We_jet1_dR;
    Histo1DPtr _h_jets_m_12;
    Histo1DPtr _h_HT;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(MC_WWJETS);

}
