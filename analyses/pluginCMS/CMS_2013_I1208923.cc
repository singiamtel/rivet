// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  // This analysis is a derived from the class Analysis:
  class CMS_2013_I1208923 : public Analysis {

  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2013_I1208923);

    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;
      declare(fs, "FS");

      // Initialize the projections
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "Jets");

      // Book histograms
      book(_h_sigma, {0., 0.5, 1., 1.5, 2., 2.5},
                     {"d01-x01-y01", "d01-x01-y02", "d01-x01-y03", "d01-x01-y04", "d01-x01-y05"});
      book(_h_invMass, {0., 0.5, 1., 1.5, 2., 2.5},
                       {"d02-x01-y01", "d02-x01-y02", "d02-x01-y03", "d02-x01-y04", "d02-x01-y05"});

    }

    // Analysis
    void analyze(const Event &event) {

      const FastJets &fJets = apply<FastJets>(event, "Jets");

      // Fill the jet pT spectra
      const Jets& jets = fJets.jetsByPt(Cuts::pt>100.*GeV && Cuts::absrap <2.5);
      for (const Jet& j : jets) {
        _h_sigma->fill(fabs(j.rapidity()), j.pT() / GeV);
      }

      // Require two jets
      const Jets& dijets = fJets.jetsByPt(Cuts::pt>30.*GeV && Cuts::absrap < 2.5);
      if (dijets.size() > 1) {
        if (dijets[0].momentum().pT() / GeV > 60.) {
          // Fill the invariant mass histogram
          double ymax = max(dijets[0].momentum().absrapidity(), dijets[1].momentum().absrapidity());
          double invMass = FourMomentum(dijets[0].momentum() + dijets[1].momentum()).mass();
          _h_invMass->fill(fabs(ymax), invMass);
        }
      }

    }


    // Scale histograms by the production cross section
    void finalize() {
      scale(_h_sigma,   crossSection() / sumOfWeights() / 2.0);
      scale(_h_invMass, crossSection() / sumOfWeights() / 2.0);
      divByGroupWidth({_h_sigma, _h_invMass});
    }

  private:

    Histo1DGroupPtr _h_sigma;
    Histo1DGroupPtr _h_invMass;

  };

  // This global object acts as a hook for the plugin system.
  RIVET_DECLARE_PLUGIN(CMS_2013_I1208923);
}
