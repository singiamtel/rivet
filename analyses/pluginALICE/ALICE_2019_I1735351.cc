// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  J/psi at 5.02 TeV
  class ALICE_2019_I1735351 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ALICE_2019_I1735351);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::pid==443), "UFS");
      book(_h_JPsi_pT,2,1,1);
      _axis = YODA::Axis<double>({ 0., 1., 2., 3., 4., 5., 7., 10. });
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_JPsi_pT->xEdges();
      // loop over J/Psi
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles()) {
        // rapidity cut
        if (p.absrap()>0.9) continue;
        _h_JPsi_pT->fill(discEdge(p.perp()));
      }
    }

    string discEdge(double val) const {
      string edge("OTHER");
      size_t idx = _axis.index(val);
      if (0 < idx && idx <= _axis.numBins())  return _edges[idx-1];
      return edge;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // factor 1/2 due folding +/- rap and 1.8 due y range
      double fact = 0.5/1.8*crossSection()/nanobarn/sumOfWeights();
      scale(_h_JPsi_pT,fact);
      for(auto & b : _h_JPsi_pT->bins()) {
        b.scaleW(1./_axis.width(b.index()));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_JPsi_pT;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ALICE_2019_I1735351);

}
