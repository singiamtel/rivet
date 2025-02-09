// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"

namespace Rivet {


  /// @brief DIS jet production cross sections at HERA
  class H1_2016_I1496981 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2016_I1496981);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(DISKinematics(), "Kinematics");

      // The final-state particles are clustered in Breit frame
      // using FastJet with the kT algorithm and a jet-radius parameter of 1.
      const DISFinalState DISfs(DISFrame::BREIT);
      FastJets jets(DISfs, fastjet::JetAlgorithm::kt_algorithm,
        fastjet::RecombinationScheme::BIpt_scheme, 1.0);
      declare(jets, "Jets");

      // Book histograms.

      // Inclusive jet pT's in Q^2 ranges.
      book(_h_ptjet_q2[0], 1, 1, 1);
      book(_h_ptjet_q2[1], 2, 1, 1);
      book(_h_ptjet_q2[2], 3, 1, 1);
      book(_h_ptjet_q2[3], 4, 1, 1);
      book(_h_ptjet_q2[4], 5, 1, 1);
      book(_h_ptjet_q2[5], 6, 1, 1);
      book(_h_ptjet_q2[6], 7, 1, 1);
      book(_h_ptjet_q2[7], 8, 1, 1);
      // High-Q^2 region.
      book(_h_high_q2, 51, 1, 1);

      // Dijet pT's.
      book(_h_ptdijet_q2[0], 9, 1, 1);
      book(_h_ptdijet_q2[1], 10, 1, 1);
      book(_h_ptdijet_q2[2], 11, 1, 1);
      book(_h_ptdijet_q2[3], 12, 1, 1);
      book(_h_ptdijet_q2[4], 13, 1, 1);
      book(_h_ptdijet_q2[5], 14, 1, 1);
      book(_h_ptdijet_q2[6], 15, 1, 1);
      book(_h_ptdijet_q2[7], 16, 1, 1);

      // Trijet pT's.
      book(_h_pttrijet_q2[0], 17, 1, 1);
      book(_h_pttrijet_q2[1], 18, 1, 1);
      book(_h_pttrijet_q2[2], 19, 1, 1);
      book(_h_pttrijet_q2[3], 20, 1, 1);
      book(_h_pttrijet_q2[4], 21, 1, 1);
      book(_h_pttrijet_q2[5], 22, 1, 1);
      book(_h_pttrijet_q2[6], 23, 1, 1);
      book(_h_pttrijet_q2[7], 24, 1, 1);
    }

    // Perform the per-event analysis
    void analyze(const Event& event) {

      // Lorentz invariant DIS quantities
      DISKinematics dis = apply<DISKinematics>(event, "Kinematics");
      double Q2 = dis.Q2();
      double y  = dis.y();
      bool isLowQ2 = (inRange(y, 0.2, 0.6) && inRange(Q2, 5.5*GeV2, 80*GeV2));
      bool isHighQ2 = (inRange(y, 0.2, 0.7) && inRange(Q2, 150.*GeV2, 15000.*GeV2));

      // Lorentz boosts for Breit and lab frames.
      const LorentzTransform breitboost = dis.boostBreit();
      const LorentzTransform labboost = breitboost.inverse();

      // Retrieve clustered jets in Breit frame, sorted by pT.
      Jets jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 4.*GeV);

      // Boost jets to lab frame.
      for (size_t i = 0; i < jets.size(); ++i) {
        jets[i].transformBy(labboost);
      }

      // Cut on Pseudorapidity in lab frame.
      // orientation 1 if hadron in "conventional" +z direction, -1 if in -z.
      const int orientation = dis.orientation();
      Jets cutJets;
      for (size_t i = 0; i < jets.size(); ++i) {
        double etaJet = jets[i].eta()*orientation;
        if ( inRange(etaJet, -1., 2.5) ) {
          cutJets += jets[i];
        }
      }

      // Boost jets back to Breit frame.
      for (size_t i = 0; i < cutJets.size(); ++i) {
        cutJets[i].transformBy(breitboost);
      }

      // Number of jets passing main cuts.
      const unsigned int N = cutJets.size();
      if (N < 1) vetoEvent;

      // Fill histograms.
      // Inclusive jets, add all jets within cuts.
      for (size_t i = 0; i < N; ++i) {

        // Doubly differential cross sections d^2 \sigma / dQ^2 dpT,
        // so need to scale filled pT's by inverse of both Q2 and pT bin widths.
        // pT is done automatically in Rivet, but Q2 must be done by
        // hand in finalize().
        const double pTJet = cutJets[i].pT();

        if (isLowQ2) {
          // Check Q2 bin
          for (size_t j = 0; j < Q2range.size()-1; ++j) {
            if ( inRange(Q2, Q2range[j], Q2range[j+1]) ) {
              _h_ptjet_q2[j]->fill(pTJet);
            }
          }
        }

        // High-Q^2 region for a pT bin.
        if (isHighQ2) {
          if (inRange(pTJet, 5.*GeV, 7.*GeV)) {
            _h_high_q2->fill(Q2);
          }
        }
      }

      // Dijets and trijets, use 2 or 3 highest pT jets.
      if (isLowQ2) {
        // Dijet <pT>
        if ( N > 1 ) {
          const double pTJet = (cutJets[0].pT() + cutJets[1].pT()) / 2.;
          for (size_t j = 0; j < Q2range.size()-1; ++j) {
            if ( inRange(Q2, Q2range[j], Q2range[j+1]) ) {
              _h_ptdijet_q2[j]->fill(pTJet);
            }
          }
        }

        // Trijet <pT>
        if ( N > 2 ) {
          const double pTJet = (cutJets[0].pT() + cutJets[1].pT() + cutJets[2].pT()) / 3.;
          for (size_t j = 0; j < Q2range.size()-1; ++j) {
            if ( inRange(Q2, Q2range[j], Q2range[j+1]) ) {
              _h_pttrijet_q2[j]->fill(pTJet);
            }
          }
        }
      }
    }

    // Normalise histograms after the run.
    // Cross sections by pT and Q2 bin widths.
    void finalize() {

      // Scaling factor to get cross section in pb.
      const double sf = crossSection()/picobarn/sumW();

      // Scale histograms by Q2 bin width.
      for (size_t i = 0; i < Q2range.size()-1; ++i) {
        scale(_h_ptjet_q2[i], sf/(Q2range[i+1]-Q2range[i]));
        scale(_h_ptdijet_q2[i], sf/(Q2range[i+1]-Q2range[i]));
        scale(_h_pttrijet_q2[i], sf/(Q2range[i+1]-Q2range[i]));
      }
      // Scale by pT bin width of 2 GeV.
      scale(_h_high_q2, sf/2.*GeV);

    }

    /// @}

  private:

    // Q2 binning required also for normalization.
    const vector<double> Q2range = { 5.5*GeV2, 8.*GeV2, 11.*GeV2, 16.*GeV2,
      22.*GeV2, 30.*GeV2, 42.*GeV2, 60.*GeV2, 80.*GeV2 };

    /// @name Histograms
    /// @{
    Histo1DPtr _h_ptjet_q2[8];
    Histo1DPtr _h_ptdijet_q2[8];
    Histo1DPtr _h_pttrijet_q2[8];
    Histo1DPtr _h_high_q2;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(H1_2016_I1496981);

}
