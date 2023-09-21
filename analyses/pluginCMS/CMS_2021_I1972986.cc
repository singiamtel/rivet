// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_2021_I1972986 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2021_I1972986);


    /// Book histograms and initialize projections:
    void init() {

      // Initialize the projections
      const FinalState fs;
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, FastJets::ANTIKT, 0.7), "JetsAK7");

      // Book sets of histograms, binned in absolute rapidity
      book(_hist_sigmaAK4, {0., 0.5, 1.0, 1.5, 2.0},
                           {"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01"});
      book(_hist_sigmaAK7, {0., 0.5, 1.0, 1.5, 2.0},
                           {"d21-x01-y01", "d22-x01-y01", "d23-x01-y01", "d24-x01-y01"});
    }


    /// Per-event analysis
    void analyze(const Event &event) {

      // AK4 jets
      const FastJets& fjAK4 = apply<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jets(Cuts::ptIn(97*GeV, 3103*GeV) && Cuts::absrap < 2.0);
      for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4->fill(j.absrap(), j.pT());
      }

      // AK7 jets
      const FastJets& fjAK7 = apply<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jets(Cuts::ptIn(97*GeV, 3103*GeV) && Cuts::absrap < 2.0);
      for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7->fill(j.absrap(), j.pT());
      }

    }


    // Finalize
    void finalize() {
      /// @note Extra division factor is the *signed* dy, i.e. 2*d|y|
      scale(_hist_sigmaAK4, 0.5*crossSection()/picobarn/sumOfWeights());
      scale(_hist_sigmaAK7, 0.5*crossSection()/picobarn/sumOfWeights());
      divByGroupWidth({_hist_sigmaAK4, _hist_sigmaAK7});
    }


    /// @name Histograms
    //@{
    Histo1DGroupPtr _hist_sigmaAK4;
    Histo1DGroupPtr _hist_sigmaAK7;
    //@}

  };


  // This global object acts as a hook for the plugin system.
  RIVET_DECLARE_PLUGIN(CMS_2021_I1972986);

}
