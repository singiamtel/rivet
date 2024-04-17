// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief CDF two-jet triply-differential cross-section
  class CDF_2001_I538041 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(CDF_2001_I538041);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      FinalState fs(Cuts::abseta < 4.2);
      declare(FastJets(fs, JetAlg::CDFJETCLU, 0.7), "Jets");

      book(_h_ET, {0.1, 0.7, 1.4, 2.1, 3.0},
                  {"d01-x01-y01", "d02-x01-y01", "d03-x01-y01", "d04-x01-y01"});
      _etaxes.resize(4);
      _etaxes[0] = YODA::Axis<double>{41.0, 47.0, 54.45, 67.2,  84.9, 100.45, 113.05, 126.0, 141.65, 162.6, 191.9, 236.75, 291.15, 350.3, 414.5};
      _etaxes[1] = YODA::Axis<double>{39.9, 46.3, 54.0,  66.85, 84.6, 100.1, 112.65, 125.55, 141.0, 161.6, 190.05, 233.7, 287.1, 343.8, 404.0};
      _etaxes[2] = YODA::Axis<double>{38.7, 45.5, 53.45, 66.15, 83.65, 98.95, 111.2, 123.75, 138.7, 158.6, 185.3, 226.75, 278.65};
      _etaxes[3] = YODA::Axis<double>{37.6, 44.2, 51.85, 63.95, 81.05, 96.1, 108.0, 120.1, 134.25, 153.1, 177.55, 205.25};
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        _edges.resize(_h_ET->numBins());
        for (const auto& b : _h_ET->bins()) {
          _edges[b.index()-1] = b->xEdges();
        }
      }

      Jets jets = apply<FastJets>(event, "Jets").jets(Cuts::Et > 10*GeV, cmpMomByEt);
      if (jets.size() < 2) vetoEvent;
      FourMomentum jet1 = jets[0].momentum();
      FourMomentum jet2 = jets[1].momentum();
      double eta1 = jet1.abseta();
      double eta2 = jet2.abseta();
      double ET1 = jet1.Et();
      double ET2 = jet2.Et();
      if (!inRange(eta1, 0.1, 0.7) || ET1 < 40.0*GeV) vetoEvent;
      if (!inRange(eta2, 0.1, 3.0)) vetoEvent;
      discfill(eta2, ET1);
      if (eta2<0.7 && ET2>40.0*GeV) discfill(eta1, ET2);
    }


    void discfill(const double eta, const double ET) {
      string edge = "OTHER";
      const size_t eta_idx = _h_ET->binAt(eta).index();
      if (eta_idx && eta_idx <= _h_ET->numBins()) {
        // eta bin is in visible range
        const size_t et_idx = _etaxes[eta_idx-1].index(ET);
        if (et_idx && et_idx <= _edges[eta_idx-1].size()) {
          // ET bin along this eta bin is also in visible range
          edge = _edges[eta_idx-1][et_idx-1];
        }
      }
      _h_ET->fill(eta, edge);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double deta1 = 1.2;
      scale(_h_ET, crossSection()/nanobarn/sumOfWeights()/deta1 / 2.0);
      divByGroupWidth(_h_ET);
      for (auto& h : _h_ET->bins()) {
        for (auto& b : h->bins()) {
          b.scaleW(1.0/_etaxes[h.index()-1].width(b.index()));
        }
      }
    }

    /// @}


  private:

    /// Histograms
    HistoGroupPtr<double,string> _h_ET;
    vector<YODA::Axis<double>> _etaxes;
    vector<vector<string>> _edges;

  };



  RIVET_DECLARE_ALIASED_PLUGIN(CDF_2001_I538041, CDF_2001_S4517016);

}
