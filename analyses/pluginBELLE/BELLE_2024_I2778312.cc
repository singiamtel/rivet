// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> D rho
  class BELLE_2024_I2778312 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2024_I2778312);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==-521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(421);
      BP.addStable(111);
      declare(BP, "BP");
      for(unsigned int ix=0;ix<2;++ix)
        book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
      	if ( ! BP.modeMatches(ix,3,mode)) continue;
        const Particle& pim = BP.decayProducts()[ix].at(-211)[0];
        const Particle& pi0 = BP.decayProducts()[ix].at( 111)[0];
        // compute the helicity angles
        LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());

        // pi pi mass and angle
        FourMomentum ppipi = pim.mom()+pi0.mom();
        ppipi = boost1.transform(ppipi);
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
        FourMomentum pPi  = boost1.transform(pim.mom());
        FourMomentum pPi3 = boost2.transform(pPi);
        double cTheta = -pPi3.p3().unit().dot(ppipi.p3().unit());
        if(cTheta>0.7) continue;
        _h[0]->fill(cTheta);
        _h[1]->fill(ppipi.mass()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
        normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode= { {-211,1}, {111,1}, {421,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2024_I2778312);

}
