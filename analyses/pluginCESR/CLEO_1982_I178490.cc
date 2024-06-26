// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"


namespace Rivet {


  /// @brief e+e- > D*+ + X
  class CLEO_1982_I178490 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1982_I178490);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::abspid==413), "UFS");
      // histos
      book(_h_sigma,1,1,1);
      book(_h_spect,2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_edges.empty()) _edges = _h_sigma->xEdges();
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        double xE = 2.*p.E()/sqrtS();
        _h_spect->fill(xE);
        for (unsigned int ix=0; ix<2; ++ix) {
          if (xE>0.5) _h_sigma->fill(_edges[0]);
          if (xE>0.7) _h_sigma->fill(_edges[1]);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_sigma,             crossSection()/sumOfWeights()/nanobarn);
      scale(_h_spect,sqr(sqrtS())*crossSection()/sumOfWeights()/microbarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_spect;
    BinnedHistoPtr<string> _h_sigma;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1982_I178490);

}
