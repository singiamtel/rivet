// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF Run I inclusive jet cross-section
  class CDF_2001_I552797 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2001_I552797);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, JetAlg::CDFJETCLU, 0.7), "Jets");
      book(_h_ET, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_ET->xEdges();
      Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 40*GeV && Cuts::abseta >= 0.1 && Cuts::abseta <= 0.7, cmpMomByEt);
      for (const Jet& jet : jets) {
        _h_ET->fill(map2string(jet.Et()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double deta = 1.2;
      scale(_h_ET, crossSection()/sumOfWeights()/deta/nanobarn);
      for (auto& b : _h_ET->bins()) {
        b.scaleW(1.0/_axis.width(b.index()));
      }
    }

    string map2string(const double value) const {
      const size_t idx = _axis.index(value);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }

    /// @}


  private:

    /// Histogram
    BinnedHistoPtr<string> _h_ET;
    vector<string> _edges;
    YODA::Axis<double> _axis{40.3, 46.3, 52.25, 58.1, 63.85, 69.5, 75.1, 80.7, 86.25, 91.75, 97.25, 102.75, 108.2,
                             113.6, 119.0, 124.4, 129.8, 135.2, 141.8, 151.05, 161.8, 172.55, 183.3, 194.1, 204.85,
                             217.8, 236.25, 257.95, 279.65, 301.3, 322.85, 347.9, 387.55, 438.25};

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2001_I552797, CDF_2001_S4563131);

}
