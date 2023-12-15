// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class CUSB_1991_I325661 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CUSB_1991_I325661);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_Bstar, "/TMP/sigma_Bstar");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      unsigned int nBstar(0);
      // Get Bottom hadrons
      const Particles bhads = ufs.particles(Cuts::abspid==513 or Cuts::abspid==523);
      // find the Bstars
      for (const Particle& p : bhads) {
        if(abs(p.pid())==513 || abs(p.pid())==523) {
          if(!p.hasDescendantWith(Cuts::pid == p.pid())) ++nBstar;
        }
      }
      if(nBstar!=0)
        _c_Bstar->fill(nBstar);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      double sig = _c_Bstar->val()*fact;
      double err = _c_Bstar->err()*fact;
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sig, err);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_Bstar;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CUSB_1991_I325661);


}
