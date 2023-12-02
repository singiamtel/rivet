// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/DileptonFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  
  /// @brief MC validation analysis for Higgs-pair events (stable Higgses)
  class MC_HHJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_HHJETS()
      : MC_JetAnalysis("MC_HHJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    /// @{

    /// Book histograms
    void init() {
      IdentifiedFinalState ifs(Cuts::abseta < 10.0 && Cuts::pT > 0*GeV);
      ifs.acceptId(25);
      declare(ifs,"IFS");

      VetoedFinalState vfs;
      vfs.addVetoPairId(25);

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

      FastJets jetpro(vfs, clusterAlgo, R);
      declare(jetpro, "Jets");

      book(_h_HH_mass ,"HH_mass", 250, 240, 4000.0);
      book(_h_HH_dR ,"HH_dR", 25, 0.5, 10.0);
      book(_h_HH_dPhi ,"HH_dPhi", 64, 0, 3.2);
      book(_h_HH_deta,"HH_deta", 50, -5, 5);
      book(_h_H_pT ,"H_pT", 50, 0, 2000.0);
      book(_h_HH_pT ,"HH_pT", 200, 0, 2000.0);
      book(_h_H_pT1 ,"H_pT1", 200, 0, 2000.0);
      book(_h_H_pT2 ,"H_pT2", 200, 0, 2000.0);
      book(_h_H_eta ,"H_eta", 50, -5.0, 5.0);
      book(_h_H_eta1 ,"H_eta1", 50, -5.0, 5.0);
      book(_h_H_eta2 ,"H_eta2", 50, -5.0, 5.0);
      book(_h_H_phi ,"H_phi", 25, 0.0, TWOPI);
      book(_h_H_jet1_deta ,"H_jet1_deta", 50, -5.0, 5.0);
      book(_h_H_jet1_dR ,"H_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {

      const IdentifiedFinalState& ifs = apply<IdentifiedFinalState>(e, "IFS");
      Particles allp = ifs.particlesByPt();
      if (allp.empty()) vetoEvent;


      FourMomentum hmom = allp[0].momentum();
      if (allp.size() > 1) {
        FourMomentum hmom2(allp[1].momentum());
        _h_HH_dR->fill(deltaR(hmom, hmom2));
        _h_HH_dPhi->fill(deltaPhi(hmom, hmom2));
        _h_HH_deta->fill(hmom.eta()-hmom2.eta());
        _h_HH_pT->fill((hmom+hmom2).pT());
        _h_HH_mass->fill((hmom+hmom2).mass());

        if (hmom.pT() > hmom2.pT()) {
          _h_H_pT1->fill(hmom.pT());
          _h_H_eta1->fill(hmom.eta());
          _h_H_pT2->fill(hmom2.pT());
          _h_H_eta2->fill(hmom2.eta());
        } else {
          _h_H_pT1->fill(hmom2.pT());
          _h_H_eta1->fill(hmom2.eta());
          _h_H_pT2->fill(hmom.pT());
          _h_H_eta2->fill(hmom.eta());
        }
      }
      _h_H_pT->fill(hmom.pT());
      _h_H_eta->fill(hmom.eta());
      _h_H_phi->fill(hmom.azimuthalAngle());


      // Get the jet candidates
      Jets jets = apply<FastJets>(e, "Jets").jetsByPt(Cuts::pT > 20*GeV);
      if (!jets.empty()) {
        _h_H_jet1_deta->fill(deltaEta(hmom, jets[0]));
        _h_H_jet1_dR->fill(deltaR(hmom, jets[0]));
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_HH_mass, crossSection()/sumOfWeights());
      scale(_h_HH_dR, crossSection()/sumOfWeights());
      scale(_h_HH_deta, crossSection()/sumOfWeights());
      scale(_h_HH_dPhi, crossSection()/sumOfWeights());
      scale(_h_H_pT, crossSection()/sumOfWeights());
      scale(_h_H_pT1, crossSection()/sumOfWeights());
      scale(_h_H_pT2, crossSection()/sumOfWeights());
      scale(_h_HH_pT, crossSection()/sumOfWeights());
      scale(_h_H_eta, crossSection()/sumOfWeights());
      scale(_h_H_eta1, crossSection()/sumOfWeights());
      scale(_h_H_eta2, crossSection()/sumOfWeights());
      scale(_h_H_phi, crossSection()/sumOfWeights());
      scale(_h_H_jet1_deta, crossSection()/sumOfWeights());
      scale(_h_H_jet1_dR, crossSection()/sumOfWeights());

      MC_JetAnalysis::finalize();
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_HH_mass, _h_HH_pT, _h_HH_dR, _h_HH_deta, _h_HH_dPhi;
    Histo1DPtr _h_H_pT, _h_H_pT1, _h_H_pT2, _h_H_eta, _h_H_eta1, _h_H_eta2, _h_H_phi;
    Histo1DPtr _h_H_jet1_deta, _h_H_jet1_dR;
    /// @}

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MC_HHJETS);

}
