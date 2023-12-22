// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF inclusive jet cross-section using the \f$ k_\perp \f$ algorithm
  class CDF_2007_I743342 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2007_I743342);


    void init() {
      // Set up projections
      const FinalState fs;
      declare(FastJets(fs, JetAlg::KT, 0.5), "JetsD05");
      declare(FastJets(fs, JetAlg::KT, 0.7), "JetsD07");
      declare(FastJets(fs, JetAlg::KT, 1.0), "JetsD10");

      // Book histos
      book(_binnedHistosD07, {0., 0.1, 0.7, 1.1, 1.6, 2.1},
                             {"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01", "d05-x01-y01"});
      book(_histoD05, 6, 1, 1);
      book(_histoD10, 7, 1, 1);
    }


    void analyze(const Event& event) {

      for (const Jet& jet : apply<JetFinder>(event, "JetsD07").jets(Cuts::pT > 54*GeV)) {
        _binnedHistosD07->fill(jet.absrap(), jet.pT());
      }

      for (const Jet& jet : apply<JetFinder>(event, "JetsD05").jets(Cuts::pT > 54*GeV)) {
        if (inRange(jet.absrap(), 0.1, 0.7))  _histoD05->fill(jet.pT());
      }

      for (const Jet& jet : apply<JetFinder>(event, "JetsD10").jets(Cuts::pT > 54*GeV)) {
        if (inRange(jet.absrap(), 0.1, 0.7))  _histoD10->fill(jet.pT());
      }
    }


    // Normalise histograms to cross-section
    void finalize() {
      const double xSec = crossSectionPerEvent()/nanobarn;

      scale(_histoD05, xSec);
      scale(_histoD10, xSec);
      // scale to xSec/yBinWidth and take into account the double yBinWidth due
      // to the absolute value of y
      scale(_binnedHistosD07, xSec/2.0);
      divByGroupWidth(_binnedHistosD07);
    }


  private:

    Histo1DGroupPtr _binnedHistosD07;

    /// Single histogram for the \f$R=0.5\f$ \f$k_\perp\f$ jets
    Histo1DPtr _histoD05;

    /// Single histogram for the \f$R=1.0\f$ \f$k_\perp\f$ jets
    Histo1DPtr _histoD10;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2007_I743342, CDF_2007_S7057202);

}
