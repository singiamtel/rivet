// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// Inclusive jet cross sections using early 13 TeV data
  /// @author Stefan von Buddenbrock <stef.von.b@cern.ch>
  class ATLAS_2016_CONF_2016_092 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_092);


    /// Bookings
    void init() {

      // Define the jets
      FastJets antiKT04Jets(FinalState(Cuts::open()),  FastJets::ANTIKT, 0.4);
      antiKT04Jets.useInvisibles();
      declare(antiKT04Jets, "antiKT04Jets");

      // Book histograms
      const vector<double> y_bins{0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0};
      book(_h_pT, y_bins);
      for (auto& b : _h_pT->bins()) {
        book(b, b.index(), 1, 1);
      }
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const Jets& jets = apply<FastJets>(event, "antiKT04Jets").jetsByPt(Cuts::pT > 100*GeV && Cuts::absrap < 3.0);
      for (const Jet& j : jets) {
        _h_pT->fill(j.absrap(), j.pT()/GeV);
      }
    }


    /// Post-run scaling
    void finalize() {
      // Divide by 2 to only get positive rapidity values
      scale(_h_pT, 0.5*crossSection()/picobarn/sumOfWeights());
    }


    /// Histograms
    Histo1DGroupPtr _h_pT;

  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(ATLAS_2016_CONF_2016_092);

}
