// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > D*+ + X
  class ARGUS_1985_I204851 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1985_I204851);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==413), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_sigma[ix],1,1,1+ix);
      }
      book(_h_spect,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges[0].empty()) {
        _edges[0] = _h_sigma[0]->xEdges();
        _edges[1] = _h_spect   ->xEdges();
      }
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const double xp = p.p3().mod()/sqrt(0.25*sqr(sqrtS())-sqr(p.mass()));
        const size_t idx = _axis.index(xp);
        if(idx&&idx<=_edges[1].size()) _h_spect->fill(_edges[1][idx-1]);
        for (unsigned int ix=0; ix<2; ++ix) {
          _h_sigma[ix]->fill(_edges[0][1]);
          if (xp>0.5) _h_sigma[ix]->fill(_edges[0][0]);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double br[2] = {0.03947,0.0822};
      const double fact = 0.677*crossSection()/sumOfWeights()/picobarn;
      for (unsigned int ix=0; ix<2; ++ix) {
        scale(_h_sigma[ix],fact*br[ix]);
      }
      scale(_h_spect,sqr(sqrtS())*crossSection()/sumOfWeights()/microbarn);
      for(auto & b : _h_spect->bins()) {
        b.scaleW(1./_axis.width(b.index()));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_spect,_h_sigma[2];
    vector<string> _edges[2];
    YODA::Axis<double> _axis{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1985_I204851);

}
