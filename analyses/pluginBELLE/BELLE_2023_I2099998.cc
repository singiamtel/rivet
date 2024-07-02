// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> K+ K- pi+
  class BELLE_2023_I2099998 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2023_I2099998);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::K0S);
      BP.addStable(  411);
      BP.addStable( -411);
      BP.addStable(  421);
      BP.addStable( -421);
      BP.addStable(10441);
      declare(BP, "BP");
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
      	book(_h[ix], 1+ix, 1, 1);
      }
      book(_c,"TMP/nB");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
        int sign = BP.decaying()[ix].pid()>0 ? 1 : -1;
        _c->fill();
        if ((sign>0 and BP.modeMatches(ix,3,mode  )) ||
            (sign<0 and BP.modeMatches(ix,3,modeCC))) {
          // 	  // momenta
          FourMomentum Kp  = BP.decayProducts()[ix].at( 321*sign)[0].mom();
          FourMomentum Km  = BP.decayProducts()[ix].at(-321*sign)[0].mom();
          FourMomentum pip = BP.decayProducts()[ix].at( 211*sign)[0].mom();
          // masses
          double mKpi = (Km+pip).mass();
          double mKK = (Kp+Km).mass();
          _h[0]->fill( mKK);
          _h[1]->fill(mKpi);
          if (mKK>1.1) continue;
          LorentzTransform boost =
            LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. mom().betaVec());
          FourMomentum pKK = boost.transform(Kp+Km);
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKK.betaVec());
          double cTheta = boost2.transform(boost.transform(Kp)).p3().unit().dot(pKK.p3().unit());
          _h[2]->fill(cTheta);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1e7/ *_c);
      normalize(_h[2], 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    CounterPtr _c;
    const map<PdgId,unsigned int> mode   = { { 321,1}, {-321,1}, { 211,1}};
    const map<PdgId,unsigned int> modeCC = { { 321,1}, {-321,1}, {-211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2023_I2099998);

}
