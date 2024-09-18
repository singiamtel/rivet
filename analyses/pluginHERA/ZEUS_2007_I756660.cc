// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DirectFinalState.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {

  /// @brief Three- and four-jet final states in photoproduction at HERA
  class ZEUS_2007_I756660 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ZEUS_2007_I756660);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Jets
      const FinalState fs;
      declare(FastJets(fs, fastjet::JetAlgorithm::kt_algorithm,
                       fastjet::RecombinationScheme::pt2_scheme, 1.0), "Jets");

      // DIS Kinematics
      declare(DISKinematics(), "Kinematics");

      // Book Histos
      book(_h["M3j"], 1, 1, 1); // T2, F3A

      book(_h["M4j"], 2, 1, 1); // T3, F3B

      // T4, F4A,4B
      book(_h["x_gamma_obs_low_M3j"], 3, 1, 1);
      book(_h["x_gamma_obs_high_M3j"], 3, 1, 2);

      // T5, F4C,4D
      book(_h["x_gamma_obs_low_M4j"], 4, 1, 1);
      book(_h["x_gamma_obs_high_M4j"], 4, 1, 2);

      // T6, F5A,5B
      book(_h["y_low_M3j"], 5, 1, 1);
      book(_h["y_high_M3j"], 5, 1, 2);

      // T7, F5C,5D
      book(_h["y_low_M4j"], 6, 1, 1);
      book(_h["y_high_M4j"], 6, 1, 2);

      // T8, F6A,6B
      book(_h["Et1_low_M3j"], 7, 1, 1);
      book(_h["Et1_high_M3j"], 7, 1, 2);

      // T9, F6A,6B
      book(_h["Et2_low_M3j"], 8, 1, 1);
      book(_h["Et2_high_M3j"], 8, 1, 2);

      // T10, F6A,6B
      book(_h["Et3_low_M3j"], 9, 1, 1);
      book(_h["Et3_high_M3j"], 9, 1, 2);

      // T11, F7A,7B
      book(_h["Et1_low_M4j"], 10, 1, 1);
      book(_h["Et1_high_M4j"], 10, 1, 2);

      // T12, F7A,7B
      book(_h["Et2_low_M4j"], 11, 1, 1);
      book(_h["Et2_high_M4j"], 11, 1, 2);

      // T13, F7A,7B
      book(_h["Et3_low_M4j"], 12, 1, 1);
      book(_h["Et3_high_M4j"], 12, 1, 2);

      // T14, F7A,7B
      book(_h["Et4_low_M4j"], 13, 1, 1);
      book(_h["Et4_high_M4j"], 13, 1, 2);

      // T15, F8A,8B
      book(_h["Eta1_low_M3j"], 14, 1, 1);
      book(_h["Eta1_high_M3j"], 14, 1, 2);

      // T16, F8A,8B
      book(_h["Eta2_low_M3j"], 15, 1, 1);
      book(_h["Eta2_high_M3j"], 15, 1, 2);

      // T17, F8A,8B
      book(_h["Eta3_low_M3j"], 16, 1, 1);
      book(_h["Eta3_high_M3j"], 16, 1, 2);

      // T18, F9A,9B
      book(_h["Eta1_low_M4j"], 17, 1, 1);
      book(_h["Eta1_high_M4j"], 17, 1, 2);

      // T19, F9A,9B
      book(_h["Eta2_low_M4j"], 18, 1, 1);
      book(_h["Eta2_high_M4j"], 18, 1, 2);

      // T20, F9A,9B
      book(_h["Eta3_low_M4j"], 19, 1, 1);
      book(_h["Eta3_high_M4j"], 19, 1, 2);

      // T21, F9A,9B
      book(_h["Eta4_low_M4j"], 20, 1, 1);
      book(_h["Eta4_high_M4j"], 20, 1, 2);

      // T22, F10A,10B
      book(_h["Cos_phi3_low_M3j"], 21, 1, 1);
      book(_h["Cos_phi3_high_M3j"], 21, 1, 2);
    }

    /// Perform the per-event analysis
    void analyze(const Event &event) {

      // DIS Kinematics & Warnings
      const DISKinematics &kin = apply<DISKinematics>(event, "Kinematics");
      if (kin.failed()) vetoEvent;
      const int orientation = kin.orientation();

      if (kin.Q2() > 1 * GeV2) vetoEvent;
      if (!inRange(kin.y(), 0.2, 0.85)) vetoEvent;

      // Jets
      const Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 6 * GeV && Cuts::etaIn(-2.4 * orientation, 2.4 * orientation), cmpMomByEt);

      // Veto event by number of jets
      if (jets.size() < 3) vetoEvent;

      // Define jet numbers
      const Jet &j1 = jets[0];
      const Jet &j2 = jets[1];
      const Jet &j3 = jets[2];

      // Define eta (with orientations)
      const double eta1 = orientation * j1.eta(), eta2 = orientation * j2.eta(), eta3 = orientation * j3.eta();

      // Get jet momenta
      FourMomentum j1_p(j1.momentum());
      FourMomentum j2_p(j2.momentum());
      FourMomentum j3_p(j3.momentum());

      // Three jet events
      if (j3.Et() > 6 * GeV) {

        // x_obs_gamma
        const double x_gamma_3j = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2) + j3.Et() * exp(-eta3)) / (2 * kin.y() * kin.beamLepton().E());

        // Invariant Mass
        const double M3j = FourMomentum(j1_p + j2_p + j3_p).mass();

        // Three Momentum
        ThreeMomentum p3 = j1_p.p3();
        ThreeMomentum p4 = j2_p.p3();
        ThreeMomentum p5 = j3_p.p3();
        ThreeMomentum p_beam = kin.beamHadron().p3() - kin.beamLepton().p3();

        Vector3 plane1 = cross(p_beam, p3);
        Vector3 plane2 = cross(p4, p5);

        const double mag_plane1 = abs(sqrt(sqr(plane1[0]) + sqr(plane1[1]) + sqr(plane1[2])));
        const double mag_plane2 = abs(sqrt(sqr(plane2[0]) + sqr(plane2[1]) + sqr(plane2[2])));

        const double norm_cos_phi3 = mag_plane1 * mag_plane2;

        const double cos_phi3 = (dot(plane1, plane2)) / norm_cos_phi3;

        _h["M3j"]->fill(M3j / GeV);

        if (inRange(M3j, 25 * GeV, 50 * GeV)) {
          _h["x_gamma_obs_low_M3j"]->fill(x_gamma_3j);
          _h["y_low_M3j"]->fill(kin.y());
          _h["Et1_low_M3j"]->fill(j1.Et());
          _h["Et2_low_M3j"]->fill(j2.Et());
          _h["Et3_low_M3j"]->fill(j3.Et());
          _h["Eta1_low_M3j"]->fill(eta1);
          _h["Eta2_low_M3j"]->fill(eta2);
          _h["Eta3_low_M3j"]->fill(eta3);
          _h["Cos_phi3_low_M3j"]->fill(cos_phi3);
        }

        if (M3j >= 50 * GeV) {
          _h["x_gamma_obs_high_M3j"]->fill(x_gamma_3j);
          _h["y_high_M3j"]->fill(kin.y());
          _h["Et1_high_M3j"]->fill(j1.Et());
          _h["Et2_high_M3j"]->fill(j2.Et());
          _h["Et3_high_M3j"]->fill(j3.Et());
          _h["Eta1_high_M3j"]->fill(eta1);
          _h["Eta2_high_M3j"]->fill(eta2);
          _h["Eta3_high_M3j"]->fill(eta3);
          _h["Cos_phi3_high_M3j"]->fill(cos_phi3);
        }
      }

      // Four jet events
      if (jets.size() > 3) {
        const Jet &j4 = jets[3];

        if (j4.Et() > 6 * GeV) {

          // Jet4 info
          FourMomentum j4_p(j4.momentum());
          const double eta4 = orientation * j4.eta();
          if (abs(eta4) > 2.4) vetoEvent;

          // x_obs_gamma
          const double x_gamma_4j = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2) + j3.Et() * exp(-eta3) + j4.Et() * exp(-eta4)) / (2 * kin.y() * kin.beamLepton().E());

          // Invariant Mass
          double M4j = FourMomentum(j1_p + j2_p + j3_p + j4_p).mass();

          _h["M4j"]->fill(M4j / GeV);

          if (inRange(M4j, 25 * GeV, 50 * GeV)) {
            _h["x_gamma_obs_low_M4j"]->fill(x_gamma_4j);
            _h["y_low_M4j"]->fill(kin.y());
            _h["Et1_low_M4j"]->fill(j1.Et());
            _h["Et2_low_M4j"]->fill(j2.Et());
            _h["Et3_low_M4j"]->fill(j3.Et());
            _h["Et4_low_M4j"]->fill(j4.Et());
            _h["Eta1_low_M4j"]->fill(eta1);
            _h["Eta2_low_M4j"]->fill(eta2);
            _h["Eta3_low_M4j"]->fill(eta3);
            _h["Eta4_low_M4j"]->fill(eta4);
          }

          if (M4j >= 50 * GeV) {
            _h["x_gamma_obs_high_M4j"]->fill(x_gamma_4j);
            _h["y_high_M4j"]->fill(kin.y());
            _h["Et1_high_M4j"]->fill(j1.Et());
            _h["Et2_high_M4j"]->fill(j2.Et());
            _h["Et3_high_M4j"]->fill(j3.Et());
            _h["Et4_high_M4j"]->fill(j4.Et());
            _h["Eta1_high_M4j"]->fill(eta1);
            _h["Eta2_high_M4j"]->fill(eta2);
            _h["Eta3_high_M4j"]->fill(eta3);
            _h["Eta4_high_M4j"]->fill(eta4);
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, crossSection() / picobarn / sumW());
    }

    /// @}

    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;


  };

  RIVET_DECLARE_PLUGIN(ZEUS_2007_I756660);

}
