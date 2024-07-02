// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief D* +/- production at 29 GeV
  class MARKII_1982_I177606 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1982_I177606);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      //Histograms
      book(_h_spect[0],2,1,1);
      book(_h_spect[1],3,1,1);
      _axis = YODA::Axis<double>(4, 0.2, 1.0);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_spect[0]->xEdges();
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==413)) {
        const double xp = 2.*p.E()/sqrtS();
        const double beta = p.p3().mod() / p.E();
        _h_spect[0]->fill(map2string(xp), 1./beta);
        _h_spect[1]->fill(map2string(xp), 1./beta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_h_spect[ix], sqr(sqrtS())*crossSection()/microbarn/sumOfWeights());
        for(auto & b : _h_spect[ix]->bins())
          b.scaleW(1./_axis.width(b.index()));
      }
    }

    ///@}

    string map2string(const double value) const {
      const size_t idx = _axis.index(value);
      if (idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER";
    }


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_spect[2];
    vector<string> _edges;
    YODA::Axis<double> _axis;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1982_I177606);

}
