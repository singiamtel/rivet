// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> D_s X
  class CLEOII_1996_I401599 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1996_I401599);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553),"UFS");
      // histograms
      book(_h, 2,1,1);
    }

    void findDecay(const Particle& parent, Particles& Ds) {
      for (const Particle& p : parent.children()) {
        if (p.abspid()==431) {
          Ds.push_back(p);
        }
        else if (p.abspid()==PID::PIPLUS ||
          p.abspid()==PID::KPLUS  ||
          p.pid()==PID::PI0 ||
          p.pid()==PID::K0S) {
          continue;
        }
        else if(!p.children().empty()) {
          findDecay(p,Ds);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles()) {
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        for (const Particle & B:p.children()) {
          if (B.abspid()!=511 && B.abspid()!=521) continue;
          Particles Ds;
          findDecay(B,Ds);
          for (const Particle & p2 : Ds) {
            const double x = boost.transform(p2.mom()).p3().mod()/sqrt(0.25*sqr(p.mass())-sqr(p2.mass()));
            _h->fill(x);
          }
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


  RIVET_DECLARE_PLUGIN(CLEOII_1996_I401599);

}
