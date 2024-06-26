// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> rho0 rho+
  class BELLE_2003_I620180 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2003_I620180);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 113);
      BP.addStable( 213);
      BP.addStable(-213);
      declare(BP, "BP");
      // histos
      book(_p[0],1,1,1);
      book(_p[1],"TMP/nW");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h_mass[ix], 2, 1, 1+ix);
      }
      book(_h_angle, 3, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
       	int sign = BP.decaying()[ix].pid()/BP.decaying()[ix].abspid();
        if ((sign== 1 && BP.modeMatches(ix,2,mode)) || (sign==-1 && BP.modeMatches(ix,2,modeCC))) {
          const Particle& rhop = BP.decayProducts()[ix].at( sign*213)[0];
          const Particle& rho0 = BP.decayProducts()[ix].at(      113)[0];
          _h_mass[0]->fill(rho0.mass());
          _h_mass[1]->fill(rhop.mass());
          // compute the helicity angle
          if (rho0.children().size()!=2) continue;
          if (rho0.children()[0].pid()!=-rho0.children()[1].pid()) continue;
          if (rho0.children()[0].abspid()!=211) continue;
          LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].mom().betaVec());
          FourMomentum prho = boost1.transform(rho0.mom());
          FourMomentum ppi  = boost1.transform(rho0.children()[0].mom());
          LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
          ppi = boost2.transform(ppi);
          const double cTheta = ppi.p3().unit().dot(prho.p3().unit());
          _h_angle->fill(cTheta);
          _p[0]->fill(-50.*(1.-5.*sqr(cTheta)));
          _p[1]->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass, 1.0, false);
      normalize(_h_angle, 1.0, false);
      scale(_p[0], 1./ *_p[1]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_angle;
    CounterPtr _p[2];
    const map<PdgId,unsigned int> mode   = { { 213,1}, {113,1} };
    const map<PdgId,unsigned int> modeCC = { {-213,1}, {113,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2003_I620180);

}
