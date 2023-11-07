// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  // Inclusive jet pT
  class CMS_2016_I1487277 : public Analysis {
  public:

    // Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_I1487277);


    // Book histograms and initialize projections:
    void init() {
      const FinalState fs;

      // Initialize the projectors:
      declare(FastJets(fs, JetAlg::ANTIKT, 0.7),"Jets");

      // Book histograms:
      book(_hist_sigma, {0., 0.5, 1., 1.5, 2., 2.5, 3., 4.7});
      for (auto& b : _hist_sigma->bins()) {
        book(b, b.index(), 1, 1);
      }

    }

    // Analysis
    void analyze(const Event &event) {
      const FastJets &fj = apply<FastJets>(event,"Jets");
      const Jets& jets = fj.jets(Cuts::ptIn(18*GeV, 5000.0*GeV) && Cuts::absrap < 5.2);

      // Fill the relevant histograms:
      for (const Jet &j : jets) {
        _hist_sigma->fill(j.absrap(), j.pT());
      }
    }

    // Finalize
    void finalize() {
      scale(_hist_sigma, crossSection()/sumOfWeights()/2.0);
    }

  private:

    Histo1DGroupPtr _hist_sigma;

  };

  // This global object acts as a hook for the plugin system.
  RIVET_DECLARE_PLUGIN(CMS_2016_I1487277);

}
