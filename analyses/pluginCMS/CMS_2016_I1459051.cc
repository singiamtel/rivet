// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Inclusive jet pT at 13 TeV
  class CMS_2016_I1459051 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2016_I1459051);


    /// Book histograms and initialize projections:
    void init() {

      // Initialize the projections
      const FinalState fs;
      declare(FastJets(fs, JetAlg::ANTIKT, 0.4), "JetsAK4");
      declare(FastJets(fs, JetAlg::ANTIKT, 0.7), "JetsAK7");

      // Book sets of histograms, binned in absolute rapidity
      book(_hist_sigmaAK7, {0., 0.5, 1., 1.5, 2., 2.5, 3.});
      book(_hist_sigmaAK4, {0., 0.5, 1., 1.5, 2., 2.5, 3.});
      for (size_t i=1; i < _hist_sigmaAK7->numBins()+1; ++i) {
        book(_hist_sigmaAK7->bin(i), i, 1, 1);
        book(_hist_sigmaAK4->bin(i), 7+i, 1, 1);
      }
      book(_hist_sigmaAK7Forward, 7, 1, 1);
      book(_hist_sigmaAK4Forward, 14, 1, 1);

    }


    /// Per-event analysis
    void analyze(const Event &event) {

      // AK4 jets
      const FastJets& fjAK4 = apply<FastJets>(event, "JetsAK4");
      const Jets& jetsAK4 = fjAK4.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      for (const Jet& j : jetsAK4) {
        _hist_sigmaAK4->fill(j.absrap(), j.pT());
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK4Forward->fill(j.pT());
      }

      // AK7 jets
      const FastJets& fjAK7 = apply<FastJets>(event, "JetsAK7");
      const Jets& jetsAK7 = fjAK7.jets(Cuts::ptIn(114*GeV, 2200.0*GeV) && Cuts::absrap < 4.7);
      for (const Jet& j : jetsAK7) {
        _hist_sigmaAK7->fill(j.absrap(), j.pT());
        if (inRange(j.absrap(), 3.2, 4.7)) _hist_sigmaAK7Forward->fill(j.pT());
      }

    }


    // Finalize
    void finalize() {
      /// @note Extra division factor is the *signed* dy, i.e. 2*d|y|
      scale(_hist_sigmaAK4, crossSection()/picobarn/sumOfWeights()/2.0);
      scale(_hist_sigmaAK7, crossSection()/picobarn/sumOfWeights()/2.0);
      scale(_hist_sigmaAK4Forward,crossSection()/picobarn/sumOfWeights()/3.0);
      scale(_hist_sigmaAK7Forward,crossSection()/picobarn/sumOfWeights()/3.0);
      divByGroupWidth({_hist_sigmaAK4, _hist_sigmaAK7});
    }


    /// @name Histograms
    /// @{
    Histo1DGroupPtr _hist_sigmaAK4;
    Histo1DGroupPtr _hist_sigmaAK7;
    Histo1DPtr _hist_sigmaAK4Forward;
    Histo1DPtr _hist_sigmaAK7Forward;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CMS_2016_I1459051);

}
