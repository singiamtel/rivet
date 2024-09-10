// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Tools/HistoGroup.hh"

namespace Rivet {

  /// @brief Energy correlators inside jets at 13 TeV
  class CMS_2024_I2760466 : public Analysis {
    public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2024_I2760466);

    /// Book histograms and initialise projections before the run
    void init() {
      // The final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      declare(FastJets(FinalState(), JetAlg::ANTIKT, 0.4), "jets");

      // Book histograms
      const vector<double> ptbins{97., 220., 330., 468., 638., 846., 1101., 1410., 1784.};
      book(_e2c, ptbins);
      book(_e3c, ptbins);
      for (unsigned int i = 0; i < ptbins.size() - 1; ++i) { // i=[0, 7]
        book(_e2c->bin(i+1), i+1, 1, 1);
        book(_e3c->bin(i+1), i+9, 1, 1);
        book(_ratios["ratio_"s+to_string(i+1)], 26 + i, 1, 1); // correspond to hepdata t[26, 33]
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Retrieve clustered jets, sorted by pT, with a minimum pT cut
      const Jets jets = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 30 * GeV);

      if (jets.size() < 2)  vetoEvent;  // at least 2 jets
      // pT and eta requirement
      for (size_t i = 0; i < 2; ++i) {
        if (jets[i].abseta() > 2.1) vetoEvent;
      }
      // back to back, diffs between leading and subleading
      const double dphi = deltaPhi(jets[0], jets[1]);
      const double dpt = jets[0].pT() - jets[1].pT();
      const double spt = jets[0].pT() + jets[1].pT();
      if (abs(dphi) <= 2 || abs(dpt) / (spt) >= 0.3) vetoEvent;

      // Calculate first-two-leading-jets e2c and e3c
      for (size_t iJet = 0; iJet < 2; ++iJet) {
        const Jet& jet = jets[iJet];
        // 1GeV cut to sub-particles inside jet
        const Particles particles = select(jet.particles(), Cuts::pT >= 1.0*GeV);
        const double E = sum(particles, Kin::E, 0.0);

        // Fill e2c and e3c histograms
        for (const Particle& p1 : particles) {
          for (const Particle& p2 : particles) {
            const double dr1 = deltaR(p1, p2);
            _e2c->fill(jet.pT(), dr1, p1.E() * p2.E() / sqr(E));
            for (const Particle& p3 : particles) {
              const double dr2 = deltaR(p3, p2);
              const double dr3 = deltaR(p3, p1);
              const double rl = max(max(dr1, dr2), dr3);
              _e3c->fill(jet.pT(), rl, p1.E() * p2.E() * p3.E() / intpow(E, 3));
            }
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_e2c);
      normalize(_e3c);
      for (const auto& b2 : _e2c->bins()) {
        if (b2->sumW() <= 0.)  continue;
        divide(_e3c->bin(b2.index()), b2, _ratios["ratio_"s+to_string(b2.index())]);
      }
    }

    private:

    Histo1DGroupPtr _e2c, _e3c;
    map<string, Estimate1DPtr> _ratios;

  };

  RIVET_DECLARE_PLUGIN(CMS_2024_I2760466);

}  // namespace Rivet
