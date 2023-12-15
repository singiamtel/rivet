// -*- C++ -*-
#include "Rivet/Analyses/MC_JETS_BASE.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z[ee]Z[mumu] + jets events
  class MC_ZZJETS : public MC_JETS_BASE {
  public:

    /// Default constructor
    MC_ZZJETS()
      : MC_JETS_BASE("MC_ZZJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      // set FS cuts from input options
      const double etaecut = getOption<double>("ABSETAEMAX", 3.5);
      const double ptecut = getOption<double>("PTEMIN", 25.);
      Cut cute = Cuts::abseta < etaecut && Cuts::pT > ptecut*GeV;
      DileptonFinder zeefinder(91.2*GeV, 0.2, cute && Cuts::abspid == PID::ELECTRON, Cuts::massIn(65*GeV, 115*GeV));
      declare(zeefinder, "ZeeFinder");
      VetoedFinalState zmminput;
      zmminput.addVetoOnThisFinalState(zeefinder);

      // set FS cuts from input options
      const double etamucut = getOption<double>("ABSETAMUMAX", 3.5);
      const double ptmucut = getOption<double>("PTMUMIN", 25.);
      Cut cutmu = Cuts::abseta < etamucut && Cuts::pT > ptmucut*GeV;
      DileptonFinder zmmfinder(zmminput, 91.2*GeV, 0.2, cutmu && Cuts::abspid == PID::MUON,
                               Cuts::massIn(65*GeV, 115*GeV));
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(zeefinder)
          .addVetoOnThisFinalState(zmmfinder);

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

      // Correlations with jets
      book(_h_ZZ_jet1_deta ,"ZZ_jet1_deta", 70, -7.0, 7.0);
      book(_h_ZZ_jet1_dR ,"ZZ_jet1_dR", 25, 1.5, 7.0);
      book(_h_Ze_jet1_dR ,"Ze_jet1_dR", 25, 0.0, 7.0);

      // Global stuff
      book(_h_HT ,"HT", logspace(100, 100.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));

      MC_JETS_BASE::init();
    }



    /// Do the analysis
    void analyze(const Event& e) {

      const DileptonFinder& zeefinder = apply<DileptonFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size() != 1) vetoEvent;

      const DileptonFinder& zmmfinder = apply<DileptonFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size() != 1) vetoEvent;

      // Z momenta
      const FourMomentum& zee = zeefinder.bosons()[0].momentum();
      const FourMomentum& zmm = zmmfinder.bosons()[0].momentum();
      const FourMomentum zz = zee + zmm;
      // Lepton momenta
      const FourMomentum& ep = zeefinder.constituents()[0].momentum();
      const FourMomentum& em = zeefinder.constituents()[1].momentum();
      const FourMomentum& mp = zmmfinder.constituents()[0].momentum();
      const FourMomentum& mm = zmmfinder.constituents()[1].momentum();

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(Cuts::pT > _jetptcut);
      if (jets.size() > 0) {
        const FourMomentum j0 = jets[0].momentum();
        _h_ZZ_jet1_deta->fill(zz.eta()-j0.eta());
        _h_ZZ_jet1_dR->fill(deltaR(zz, j0));
        _h_Ze_jet1_dR->fill(deltaR(ep, j0));
      }

      const double HT = sum(jets, Kin::pT, ep.pT() + em.pT() + mp.pT() + mm.pT());
      if (HT > 0.0) _h_HT->fill(HT/GeV);

      MC_JETS_BASE::analyze(e);
    }


    /// Finalize
    void finalize() {
      const double s = crossSection()/picobarn/sumOfWeights();
      scale(_h_ZZ_jet1_deta, s);
      scale(_h_ZZ_jet1_dR, s);
      scale(_h_Ze_jet1_dR, s);
      scale(_h_HT, s);
      MC_JETS_BASE::finalize();
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_ZZ_jet1_deta;
    Histo1DPtr _h_ZZ_jet1_dR;
    Histo1DPtr _h_Ze_jet1_dR;
    Histo1DPtr _h_HT;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(MC_ZZJETS);

}
