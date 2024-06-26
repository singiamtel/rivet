// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> lambdabar/Sigmabar0 p pi-
  class BELLE_2023_I2663731 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2663731);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable( 3122);
      B0.addStable(-3122);
      B0.addStable( 3212);
      B0.addStable(-3212);
      declare(B0, "B0");
      // histos
      for (unsigned int ix=0;ix<2;++ix) {
        book(_h_m    [ix],1,1,1+ix);
        book(_h_theta[ix],2,1,1+ix);
      }
      book(_cB,"/TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for (unsigned int ix=0; ix<B0.decaying().size(); ++ix) {
        _cB->fill();
        unsigned int imode=0;
        if      (B0.modeMatches(ix,3,mode1)) imode=0;
        else if (B0.modeMatches(ix,3,mode2)) imode=1;
        else continue;
        const Particle & pp     = B0.decayProducts()[ix].at( 2212)[0];
        const Particle & ppim   = B0.decayProducts()[ix].at(- 211)[0];
        const Particle & LamBar = B0.decayProducts()[ix].at(imode==0 ? -3122 : -3212)[0];
        double mBB =  (pp.mom()+LamBar.mom()).mass();
        _h_m[imode]->fill(mBB);
        if (mBB>2.8) continue;
        // boost to B rest frame
        LorentzTransform boost =
          LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].mom().betaVec());
        FourMomentum pLam    = boost.transform(LamBar.mom());
        FourMomentum pProton = boost.transform(pp    .mom());
        Vector3 axis1 =  boost.transform(ppim.mom()).p3().unit();
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta((pLam+pProton).betaVec());
        Vector3 axis2 = boost2.transform(pProton).p3().unit();
        _h_theta[imode]->fill(axis1.dot(axis2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_m    , 1.e6/ *_cB);
      scale(_h_theta, 1.e6/ *_cB);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_m[2],_h_theta[2];
    CounterPtr _cB;
    const map<PdgId,unsigned int> mode1   = { { 2212,1},{-3122,1}, {-211,1}};
    const map<PdgId,unsigned int> mode2   = { { 2212,1},{-3212,1}, {-211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2663731);

}
