// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B+ -> psi(3770) K+
  class BELLE_2004_I623994 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2004_I623994);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      // histos
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      UnstableParticles ufs = apply<UnstableParticles>(event,"UFS");
      for (const Particle& B : ufs.particles()) {
      	if (B.children().size()!=2) continue;
      	int sign = B.pid()/B.abspid();
        // boost to B rest frame
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B.mom().betaVec());
        Particle resonance;
        if (B.children()[0].pid()==30443 && B.children()[1].pid()==sign*321) {
          resonance=B.children()[0];
        }
        else if (B.children()[1].pid()==30443 && B.children()[0].pid()==sign*321) {
          resonance=B.children()[1];
        }
        else {
          continue;
        }
        if (resonance.children().size()!=2) continue;
        FourMomentum pRes = boost1.transform(resonance.mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pRes.betaVec());
        Vector3 axis1 = pRes.p3().unit();
        Particle child;
        if (resonance.children()[0].pid()!=-resonance.children()[1].pid()) continue;
        if (resonance.children()[0].abspid()!=411 && resonance.children()[0].abspid()!=421) continue;
        child = resonance.children()[0].pid()>0 ? resonance.children()[0] : resonance.children()[1];
        Vector3 axis2 = boost2.transform(boost1.transform(child.mom())).p3().unit();
      	_h->fill(axis1.dot(axis2));
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


  RIVET_DECLARE_PLUGIN(BELLE_2004_I623994);

}
