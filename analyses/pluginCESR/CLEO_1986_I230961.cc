// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi spectrum in Upsilon(4S) decays
  class CLEO_1986_I230961 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1986_I230961);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553), "UFS");
      // histos
      book(_h,2,1,1);
    }

    /// Recursively walk the decay tree to find decay products of @a p
    void findDecayProducts(Particle mother, Particles& unstable) {
      for(const Particle& p: mother.children()) {
        const int id = abs(p.pid());
        if (id == 443) {
          unstable.push_back(p);
        }
        else if(!p.children().empty()) {
          findDecayProducts(p, unstable);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& ups : ufs.particles()) {
        Particles unstable;
        // Find the decay products we want
        findDecayProducts(ups,unstable);
        LorentzTransform cms_boost;
        if (ups.p3().mod() > 0.001) {
          cms_boost = LorentzTransform::mkFrameTransformFromBeta(ups.mom().betaVec());
        }
        for (const Particle& p : unstable) {
          const double modp = cms_boost.transform(p.mom()).p3().mod();
          _h->fill(modp);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1986_I230961);

}
