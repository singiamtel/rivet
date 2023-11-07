// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for W^+[enu]W^-[munu] + jets events
  class MC_WWJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_WWJETS()
      : MC_JetAnalysis("MC_WWJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      FinalState fs;

      // set FS cuts from input options
      const double etaecut = getOption<double>("ABSETAEMAX", 3.5);
      const double ptecut = getOption<double>("PTEMIN", 25.);

      Cut cute = Cuts::abseta < etaecut && Cuts::pT > ptecut*GeV;

      WFinder wenufinder(fs, cute, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wenufinder, "WenuFinder");

      VetoedFinalState wmnuinput;
      wmnuinput.addVetoOnThisFinalState(wenufinder);

      // set FS cuts from input options
      const double etamucut = getOption<double>("ABSETAMUMAX", 3.5);
      const double ptmucut = getOption<double>("PTMUMIN", 25.);

      Cut cutmu = Cuts::abseta < etamucut && Cuts::pT > ptmucut*GeV;
      
      WFinder wmnufinder(wmnuinput, cutmu, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wmnufinder, "WmnuFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(wenufinder)
          .addVetoOnThisFinalState(wmnufinder);

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

      FastJets jetpro(jetinput, clusterAlgo, R);
      declare(jetpro, "Jets");

      // correlations with jets
      book(_h_WW_jet1_deta ,"WW_jet1_deta", 70, -7.0, 7.0);
      book(_h_WW_jet1_dR ,"WW_jet1_dR", 25, 1.5, 7.0);
      book(_h_We_jet1_dR ,"We_jet1_dR", 25, 0.0, 7.0);

      // global stuff
      book(_h_HT ,"HT", logspace(100, 100.0, 0.5*(sqrtS()>0.?sqrtS():14000.)));
      book(_h_jets_m_12 ,"jets_m_12", logspace(100, 1.0, 0.25*(sqrtS()>0.?sqrtS():14000.)));

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event& e) {
      const double weight = 1.0;

      const WFinder& wenufinder = apply<WFinder>(e, "WenuFinder");
      if (wenufinder.bosons().size() !=1 ) vetoEvent;

      const WFinder& wmnufinder = apply<WFinder>(e, "WmnuFinder");
      if (wmnufinder.bosons().size() !=1 ) vetoEvent;

      FourMomentum wenu = wenufinder.bosons()[0].momentum();
      FourMomentum wmnu = wmnufinder.bosons()[0].momentum();
      FourMomentum ww = wenu + wmnu;
      // find leptons
      FourMomentum ep = wenufinder.leptons()[0].momentum();
      FourMomentum enu = wenufinder.neutrinos()[0].momentum();
      FourMomentum mm = wmnufinder.leptons()[0].momentum();
      FourMomentum mnu = wmnufinder.neutrinos()[0].momentum();

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size() > 0) {
        _h_WW_jet1_deta->fill(ww.eta()-jets[0].eta(), weight);
        _h_WW_jet1_dR->fill(deltaR(ww, jets[0].momentum()), weight);
        _h_We_jet1_dR->fill(deltaR(ep, jets[0].momentum()), weight);
      }

      double HT = ep.pT() + mm.pT() + FourMomentum(enu+mnu).pT();
      for (const Jet& jet : jets) HT += jet.pT();
      if (HT > 0.0) _h_HT->fill(HT/GeV, weight);

      if (jets.size() > 1) {
        const FourMomentum jet1 = jets[0].momentum();
        const FourMomentum jet2 = jets[1].momentum();
        _h_jets_m_12->fill((jet1+jet2).mass()/GeV, weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/picobarn/sumOfWeights();
      scale(_h_WW_jet1_deta, norm);
      scale(_h_WW_jet1_dR, norm);
      scale(_h_We_jet1_dR, norm);
      scale(_h_jets_m_12, norm);
      scale(_h_HT, norm);
      MC_JetAnalysis::finalize();
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



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_WWJETS);

}
