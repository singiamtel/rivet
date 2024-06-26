// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> rho0 pi0
  class BELLE_2006_I688850 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I688850);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 111);
      declare(B0, "B0");
      // histos
      for(unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      for(unsigned int ix=0; ix<B0.decaying().size(); ++ix) {
      	int sign = B0.decaying()[ix].pid()/B0.decaying()[ix].abspid();
      	if (!B0.modeMatches(ix,3,mode)) continue;
       	const Particle& pip = B0.decayProducts()[ix].at( sign*211)[0];
       	const Particle& pim = B0.decayProducts()[ix].at(-sign*211)[0];
       	const Particle& pi0 = B0.decayProducts()[ix].at(      111)[0];
        // require pi+- pi0 mass >1.1 GeV
        if((pip.mom()+pi0.mom()).mass()<1.1) continue;
        if((pim.mom()+pi0.mom()).mass()<1.1) continue;
        FourMomentum pRho = pip.mom()+pim.mom();
        _h[0]->fill(pRho.mass());
       	// compute the helicity angle
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].mom().betaVec());
      	pRho = boost1.transform(pRho);
      	FourMomentum ppi  = boost1.transform(pip.mom());
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pRho.betaVec());
        ppi = boost2.transform(ppi);
        const double cTheta = ppi.p3().unit().dot(pRho.p3().unit());
        _h[1]->fill(cTheta);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode = { { 211,1}, {-211,1}, { 111,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I688850);

}
