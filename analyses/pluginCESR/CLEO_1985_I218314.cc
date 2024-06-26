// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> phi X
  class CLEO_1985_I218314 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1985_I218314);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553),"UFS");
      // histograms
      book(_h,4,1,1);
      book(_c,"TMP/nB");
    }

    void findDecay(const Particle& parent, Particles& phi) {
      for (const Particle& p : parent.children()) {
        if (p.abspid()==333) {
          phi.push_back(p);
        }
        else if (p.abspid()==PID::PIPLUS ||
                 p.abspid()==PID::KPLUS  ||
                 p.pid()==PID::PI0 ||
                 p.pid()==PID::K0S)
          continue;
        else if (!p.children().empty()) {
          findDecay(p,phi);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        for (const Particle& B:p.children()) {
          if (B.abspid()!=511 && B.abspid()!=521) continue;
          _c->fill();
          Particles phi;
          findDecay(B,phi);
          for (const Particle & p2 : phi) {
            const double pmod = boost.transform(p2.mom()).p3().mod();
            _h->fill(pmod);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,1./ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1985_I218314);

}
