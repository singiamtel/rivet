// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brieftau -> hadrons
  class CLEO_1995_I396711 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1995_I396711);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(Cuts::abspid==PID::TAU), "UFS");
      book(_h, 1, 1, 1);
    }


    void findDecayProducts(const Particle& mother, Particles& hadrons) {
      for (const Particle& p : mother.children()) {
        long id = p.abspid();
        if (id == PID::PI0 || id==311 || id==310 || id==130) {
          hadrons.push_back(p);
        }
        else if (!p.children().empty()) {
          findDecayProducts(p, hadrons);
        }
        else if (id>16) {
          hadrons.push_back(p);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles(Cuts::abspid==PID::TAU)) {
        Particles hadrons;
        // find the decay products we want
        findDecayProducts(p, hadrons);
        if (hadrons.size()<=1) continue;
        FourMomentum ptotal;
        for (const Particle& p2 : hadrons) ptotal+=p2.momentum();
        _h->fill(ptotal.mass2());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1995_I396711);

}
