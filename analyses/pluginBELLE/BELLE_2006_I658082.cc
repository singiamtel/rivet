// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> Lambdac+ pbar pi-
  class BELLE_2006_I658082 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2006_I658082);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      // treat Lambda_c, Sigma_c and sigma_c* as stable
      BP.addStable( 4122);
      BP.addStable(-4122);
      declare(BP, "BP");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
      book(_h[2],2,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0;ix<BP.decaying().size();++ix) {
       	int sign = 1;
       	if       (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,mode  )) sign = 1;
       	else if  (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,modeCC)) sign =-1;
        else continue;
        const Particle& lamC = BP.decayProducts()[ix].at( sign*4122)[0];
        const Particle& pbar = BP.decayProducts()[ix].at(-sign*2212)[0];
        const Particle& pim  = BP.decayProducts()[ix].at(-sign*211 )[0];
        double mLamCPi = (lamC.mom()+pim .mom()).mass();
        double mPbarPi = (pbar.mom()+pim .mom()).mass();
        double mLamCPB = (lamC.mom()+pbar.mom()).mass();
        if (mLamCPi>2.6 && mLamCPB>3.5) _h[0]->fill(mPbarPi);
        if (mLamCPi>2.6 && mPbarPi>1.6) _h[1]->fill(mLamCPB);
        if (mLamCPB<3.6) {
          LorentzTransform boostB = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());
          FourMomentum pLamC   = boostB.transform(lamC.mom());
          FourMomentum pLamCPB = boostB.transform(lamC.mom()+pbar.mom());
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pLamCPB.betaVec());
          pLamC   = boost2.transform(pLamC);
          double cTheta = pLamC.p3().unit().dot(pLamCPB.p3().unit());
          _h[2]->fill(cTheta);
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
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode   = { { 4122,1},{-2212,1}, {-211,1}};
    const map<PdgId,unsigned int> modeCC = { {-4122,1},{ 2212,1}, { 211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2006_I658082);

}
