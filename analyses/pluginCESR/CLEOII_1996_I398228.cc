// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B -> eta X
  class CLEOII_1996_I398228 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1996_I398228);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(UnstableParticles(Cuts::pid==300553),"UFS");
      // histograms
      book(_w_ups4,"/TMP/w_ups4");
      book(_h,1,1,1);
    }

    void findDecay(const Particle& parent, Particles& eta) {
      for (const Particle& p : parent.children()) {
        if (p.pid()==PID::ETA) {
          eta.push_back(p);
        }
        else if (p.abspid()==PID::PIPLUS ||
                 p.abspid()==PID::KPLUS  ||
                 p.pid()==PID::PI0 ||
                 p.pid()==PID::K0S) {
          continue;
        }
        else if(!p.children().empty()) {
          findDecay(p, eta);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        _w_ups4->fill();
        const LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(p.mom().betaVec());
        for (const Particle& B : p.children()) {
          if (B.abspid()!=511 && B.abspid()!=521) continue;
          Particles eta;
          findDecay(B, eta);
          for (const Particle& p2 : eta) {
            double pEta = boost.transform(p2.mom()).p3().mod();
            _h->fill(pEta);
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 2 B's per upsilon and in percent
      scale(_h, 50./ *_w_ups4);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _w_ups4;
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1996_I398228);

}
