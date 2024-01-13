// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/DISKinematics.hh"
#include "Rivet/Projections/DISFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief H1 dijet photoproduction
  ///
  /// @note Cleaning cuts on event pT/sqrt(Et) and y_e are not needed in MC analysis.
  class H1_2006_I711847 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(H1_2006_I711847);

    /// @name Analysis methods
    //@{

    // Book projections and histograms
    void init() {

      /// @todo Acceptance
      const FinalState & fs = declare(DISFinalState(DISFrame::LAB), "FS");
      declare(FastJets(fs, JetAlg::KT, 1.0), "Jets");

      // Projections
      declare(DISKinematics(), "Kinematics");

      // Table 1
      book(_h_costh[0], 1, 1, 1);
      book(_h_costh[1], 1, 1, 2);
      // Table 2
      book(_h_costh_mjj[0], 2, 1, 1);
      book(_h_costh_mjj[1], 2, 1, 2);
      // Table 3
      book(_h_xgamma[0], 3, 1, 1);
      book(_h_xgamma[1], 3, 1, 2);
      // Table 4
      book(_h_xproton[0], 4, 1, 1);
      book(_h_xproton[1], 5, 1, 1);
      book(_h_xproton[2], 6, 1, 1);
      book(_h_xproton[3], 7, 1, 1);
      book(_h_xproton[4], 8, 1, 1);
      book(_h_xproton[5], 9, 1, 1);
      // Table 5
      book(_h_etjet1[0], 10, 1, 1);
      book(_h_etjet1[1], 11, 1, 1);
      book(_h_etjet1[2], 12, 1, 1);
      book(_h_etjet1[3], 13, 1, 1);
      book(_h_etjet1[4], 14, 1, 1);
      book(_h_etjet1[5], 15, 1, 1);
    }

    // Do the analysis
    void analyze(const Event& event) {

      // Determine event orientation, since coord system is for +z = proton direction
      const ParticlePair bs = event.beams();
      if (bs.first.pid() != PID::POSITRON && bs.second.pid() != PID::POSITRON) vetoEvent;
      const Particle& bpositron = (bs.first.pid() == PID::POSITRON ? bs.first : bs.second);
      if (bs.first.pid() != PID::PROTON && bs.second.pid() != PID::PROTON) vetoEvent;
      const Particle& bproton = (bs.first.pid() == PID::PROTON) ? bs.first : bs.second;
      const int orientation = sign(bproton.momentum().pz());
      MSG_DEBUG("Beam proton = " << bproton.mom() << " GeV => orientation = " << orientation);

      // Jet selection
      const Jets jets = apply<FastJets>(event, "Jets") \
        .jets(Cuts::Et > 15*GeV && Cuts::etaIn(-0.5*orientation, 2.75*orientation), cmpMomByEt);
      MSG_DEBUG("Jet multiplicity = " << jets.size());
      if (jets.size() < 2) vetoEvent;
      const Jet& j1 = jets[0];
      const Jet& j2 = jets[1];
      if (j1.Et() < 25*GeV) vetoEvent;

      // eta and cos(theta*) computation
      const double eta1 = orientation*j1.eta(), eta2 = orientation*j2.eta();
      const double etadiff = eta1 - eta2;
      const double costhetastar = tanh(0.5*etadiff);

      // DIS kinematics
      const DISKinematics& dk = apply<DISKinematics>(event, "Kinematics");
      double q2  = dk.Q2();
      double y   = dk.y();

      if (q2 > 1.) vetoEvent;

      // Computation and cut on inelasticity
      if (!inRange(y, 0.1, 0.9)) vetoEvent;

      // Computation of x_y^obs
      const double xphoton = (j1.Et() * exp(-eta1) + j2.Et() * exp(-eta2)) / (2*y*bpositron.E());
      const double xproton = (j1.Et() * exp( eta1) + j2.Et() * exp( eta2)) / (2*bproton.E());

      const size_t i_xphoton = (xphoton < 0.8) ? 0 : 1;
      const size_t i_xproton = (xproton < 0.1) ? 0 : 1;

      // Calculate the invariant mass of the dijet as in the paper
      const double mjj = sqrt( 2.*j1.Et()*j2.Et()*( cosh(j1.eta() - j2.eta()) - cos(j1.phi() - j2.phi()) ) );

      // Fill histograms
      // T1
      _h_costh[i_xphoton]->fill(abs(costhetastar));
      // T2
      if (mjj > 65*GeV) _h_costh_mjj[i_xphoton]->fill(abs(costhetastar));
      // T3
      _h_xgamma[i_xproton]->fill(xphoton);

      // T4
      if (eta1 < 1 && eta2 < 1 ) _h_xproton[3 * i_xphoton]->fill(xproton);
      if (eta1 > 1 && eta2 < 1 ) _h_xproton[1 + 3 * i_xphoton]->fill(xproton);
      if (eta2 > 1 && eta1 < 1 ) _h_xproton[1 + 3 * i_xphoton]->fill(xproton);
      if (eta1 > 1 && eta2 > 1 ) _h_xproton[2 + 3 * i_xphoton]->fill(xproton);
      // T5
      if (eta1 < 1 && eta2 < 1 ) _h_etjet1[3 * i_xphoton]->fill(j1.Et()/GeV);
      if (eta1 > 1 && eta2 < 1 ) _h_etjet1[1 + 3 * i_xphoton]->fill(j1.Et()/GeV);
      if (eta2 > 1 && eta1 < 1 ) _h_etjet1[1 + 3 * i_xphoton]->fill(j1.Et()/GeV);
      if (eta1 > 1 && eta2 > 1 ) _h_etjet1[2 + 3 * i_xphoton]->fill(j1.Et()/GeV);
    // }
    }

    // Finalize
    void finalize() {
      const double sf = crossSection()/picobarn/sumOfWeights();
      scale(_h_costh, sf);
      scale(_h_costh_mjj, sf);
      scale(_h_xgamma, sf);
      scale(_h_xproton, sf);
      scale(_h_etjet1, sf);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_costh[2], _h_costh_mjj[2], _h_xgamma[2], _h_xproton[6], _h_etjet1[6];
    //@}

  };


  RIVET_DECLARE_PLUGIN(H1_2006_I711847);

}
