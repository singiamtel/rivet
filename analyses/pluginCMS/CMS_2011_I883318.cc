// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B+ at 7 TeV
  class CMS_2011_I883318 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2011_I883318);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      declare(UnstableParticles(), "UFS");
      // histograms
      book(_h_total, 1, 1, 1);
      book(_h_pT   , 2, 1, 1);
      book(_h_y    , 3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges[0].empty()) {
        _edges[0] = _h_total->xEdges();
        _edges[1] = _h_pT->xEdges();
      }
      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // loop over B+ states
      for (const Particle& p : ufs.particles(Cuts::abspid==521)) {
        // cuts on pT and rapidity
        const double y  = p.absrap();
        const double pT = p.perp();
        if (pT<5. || y>2.4) continue;
        _h_total->fill(_edges[0][0]);
        const size_t idx = _axis.index(pT);
        const string edge = (idx&&idx<_edges[1].size()) ? _edges[1][idx-1] : "OTHER"s;
        _h_pT->fill(edge);
        _h_y->fill(y);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = 0.5*crossSection()/microbarn/sumOfWeights();
      scale(_h_total, fact);
      scale(_h_pT, fact);
      for (auto& b : _h_pT->bins()) {
        const size_t idx = b.index();
        b.scaleW(1.0/_axis.width(idx));
      }
      // 0.5 from +/- y
      scale(_h_y, 0.5*fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_total,_h_pT;
    Histo1DPtr _h_y;
    vector<string> _edges[2];
    YODA::Axis<double> _axis{5.0,10.0,13.0,17.0,24.0,30.0};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMS_2011_I883318);

}
