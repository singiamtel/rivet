// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class CLEO_1991_I29927 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1991_I29927);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_B, "/TMP/sigma_B");
      book(_c_Bstar, "/TMP/sigma_Bstar");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      unsigned int nBstar(0);
      // Get Bottom hadrons
      const Particles bhads = select(ufs.particles(), isBottomHadron);
      // find the Bstars
      for (const Particle& p : bhads) {
        if(abs(p.pid())==513 || abs(p.pid())==523) {
          if(!p.hasDescendantWith(Cuts::pid == p.pid())) ++nBstar;
        }
      }
      if(!bhads.empty())
        _c_B->fill();
      if(nBstar!=0)
        _c_Bstar->fill(nBstar);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=1;ix<3;++ix) {
        double sig(0.),err(0.);
        if(ix==1) {
          sig = _c_B->val()*fact;
          err = _c_B->err()*fact;
        }
        else {
          sig = _c_Bstar->val()*fact;
          err = _c_Bstar->err()*fact;
        }
        Estimate1DPtr mult;
        book(mult, ix, 1, 1);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS()/GeV, b.xMid(), b.xMax())) {
            b.set(sig, err);
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_B, _c_Bstar;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1991_I29927);


}
