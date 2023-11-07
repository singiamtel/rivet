// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run II inclusive jet cross-section using the Midpoint algorithm.
  ///
  /// The analysis includes 1.1fb^-1 of CDF data and is the first with a
  /// cone algorithm to include the forward region of the detector.
  /// arXiv:0807.2204 to be published in PRD
  class CDF_2008_S7828950 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2008_S7828950);


    /// @name Analysis methods
    /// @{

    // Book histos and set counters for number of events passed in each one
    void init() {
      const FinalState fs;
      declare(FastJets(fs, JetAlg::CDFMIDPOINT, 0.7), "JetsM07");

      book(_binnedHistosR07, {0., 0.1, 0.7, 1.1, 1.6, 2.1},
                             {"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01", "d05-x01-y01"});
    }


    // Do the analysis
    void analyze(const Event& event) {
      for (const Jet& jet : apply<FastJets>(event, "JetsM07").jets(Cuts::pT > 62*GeV)) {
        _binnedHistosR07->fill(jet.absrap(), jet.pT()/GeV);
      }
    }


    // Normalise histograms to cross-section
    void finalize() {
      scale(_binnedHistosR07, crossSection()/nanobarn/sumOfWeights()/2.0);
      divByGroupWidth(_binnedHistosR07);
    }

    /// @}


  private:

    /// Histograms in different eta regions
    Histo1DGroupPtr _binnedHistosR07;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2008_S7828950, CDF_2008_I790693);

}
