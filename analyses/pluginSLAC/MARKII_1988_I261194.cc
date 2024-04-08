// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta production at 29 GeV
  class MARKII_1988_I261194 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1988_I261194);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_spect,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_spect->xEdges();
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==221)) {
        const double xp = 2.*p.E()/sqrtS();
        const double beta = p.p3().mod() / p.E();
        _h_spect->fill(map2string(xp), 1./beta);
      }
    }

    string map2string(const double value) const {
      const size_t idx = _axis.index(value);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_spect, sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
        for(auto & b: _h_spect->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axis.width(idx));
        }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_spect;
    vector<string> _edges;
    YODA::Axis<double> _axis{0.2, 0.3, 0.35, 0.4, 0.45, 0.525, 0.625, 0.725};
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1988_I261194);

}
