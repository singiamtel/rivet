// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"

namespace Rivet {


  /// @brief Inclusive jet cross sections, differential in
  /// jet transverse energy E_T.
  class H1_2002_I588263 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2002_I588263);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISKinematics(), "Kinematics");

      // The final-state particles are clustered in Breit frame
      // using FastJet with the kT algorithm and a jet-radius parameter of 1.
      const DISFinalState DISfs(DISFrame::BREIT);
      FastJets jets(DISfs, JetAlg::KT, 1.0);
      declare(jets, "Jets");

      // Book histograms.
      // Transverse jet energies separated into pseudorapidity ranges.
      book(_h_etjet_eta[0], 1, 1, 1);
      book(_h_etjet_eta[1], 1, 1, 2);
      book(_h_etjet_eta[2], 1, 1, 3);

      // Transverse jet energies separated into Q2 ranges.
      book(_h_etjet_q2[0], 2, 1, 1);
      book(_h_etjet_q2[1], 2, 1, 2);
      book(_h_etjet_q2[2], 2, 1, 3);
      book(_h_etjet_q2[3], 2, 1, 4);
      book(_h_etjet_q2[4], 2, 1, 5);

      // ET^2 by Q^2, separated into pseudorapidity ranges.
      book(_h_et2q2jet_eta[0], 3, 1, 1);
      book(_h_et2q2jet_eta[1], 4, 1, 1);
      book(_h_et2q2jet_eta[2], 5, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Lorentz invariant DIS quantities
      DISKinematics dis = apply<DISKinematics>(event, "Kinematics");
      double Q2 = dis.Q2();
      double y  = dis.y();

      // Kinematic cuts on virtuality and inelasticity.
      if ( !inRange(Q2, 5.*GeV2, 100.*GeV2) ) vetoEvent;
      if ( !inRange(y, 0.2, 0.6) )            vetoEvent;

      // Lorentz boosts for Breit and lab frames.
      const LorentzTransform breitboost = dis.boostBreit();
      const LorentzTransform labboost = breitboost.inverse();

      // Retrieve clustered jets in Breit frame, sorted by pT.
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::Et > 5*GeV);

      // Boost jets to lab frame.
      for (size_t i = 0; i < jets.size(); ++i) {
        jets[i].transformBy(labboost);
      }

      // Cut on Pseurdorapidity in lab frame.
      // 1 if hadron in "conventional" +z direction, -1 if in -z.
      const int orientation = dis.orientation();
      Jets labJets;
      for (size_t i = 0; i < jets.size(); ++i) {
        double etaJet = jets[i].eta()*orientation;
        if ( inRange(etaJet, -1., 2.8) ) {
          labJets += jets[i];
        }
      }

      // Boost jets back to Breit frame.
      Jets breitJets = labJets;
      for (size_t i = 0; i < breitJets.size(); ++i){
        breitJets[i].transformBy(breitboost);
      }

      // Fill histograms.
      for (size_t i = 0; i < labJets.size(); ++i) {

        // Jets are in the same order in both frames.
        double etaJet = labJets[i].eta()*orientation;
        double eTJet  = breitJets[i].Et();
        double eT2Q2  = eTJet*eTJet / Q2;

        // ET in 3 different eta ranges
        if ( inRange(etaJet, -1., 0.5) ) {
          _h_etjet_eta[0]->fill(eTJet);
          _h_et2q2jet_eta[0]->fill(eT2Q2);
        } else if ( inRange(etaJet, 0.5, 1.5) ) {
          _h_etjet_eta[1]->fill(eTJet);
          _h_et2q2jet_eta[1]->fill(eT2Q2);
        } else if ( inRange(etaJet, 1.5, 2.8) ) {
          _h_etjet_eta[2]->fill(eTJet);
          _h_et2q2jet_eta[2]->fill(eT2Q2);
        }
        // ET of forward region, in 5 different Q2 ranges
        if ( inRange(etaJet, 1.5, 2.8) ){
          if (5*GeV2 < Q2 && Q2 < 10*GeV2) {
            _h_etjet_q2[0]->fill(eTJet);
          } else if (10*GeV2 < Q2 && Q2 < 20*GeV2) {
            _h_etjet_q2[1]->fill(eTJet);
          } else if (20*GeV2 < Q2 && Q2 < 35*GeV2) {
            _h_etjet_q2[2]->fill(eTJet);
          } else if (35*GeV2 < Q2 && Q2 < 70*GeV2) {
            _h_etjet_q2[3]->fill(eTJet);
          } else if (70*GeV2 < Q2 && Q2 < 100*GeV2) {
            _h_etjet_q2[4]->fill(eTJet);
          }
        }
      }
    }

    /// Normalise histograms after the run
    void finalize() {
      // Scaling factor
      const double sf = crossSection()/picobarn/sumW();

      scale(_h_etjet_eta[0], sf);
      scale(_h_etjet_eta[1], sf);
      scale(_h_etjet_eta[2], sf);

      scale(_h_etjet_q2[0], sf);
      scale(_h_etjet_q2[1], sf);
      scale(_h_etjet_q2[2], sf);
      scale(_h_etjet_q2[3], sf);
      scale(_h_etjet_q2[4], sf);

      scale(_h_et2q2jet_eta[0], sf);
      scale(_h_et2q2jet_eta[1], sf);
      scale(_h_et2q2jet_eta[2], sf);
    }

    /// @}

  private:

    /// @name Histograms
    /// @{
    Histo1DPtr _h_etjet_eta[3];
    Histo1DPtr _h_etjet_q2[5];
    Histo1DPtr _h_et2q2jet_eta[3];
    /// @}
  };


  RIVET_DECLARE_PLUGIN(H1_2002_I588263);

}
