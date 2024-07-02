// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"


namespace Rivet {


  /// @brief J/psi -> p K- Lambdabar
  class BESII_2004_I650570 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2004_I650570);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable( PID::PI0);
      psi.addStable( PID::K0S);
      psi.addStable( PID::LAMBDA);
      psi.addStable(-PID::LAMBDA);
      declare(psi, "psi");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0; ix<psi.decaying().size(); ++ix) {
        int sign=1;
      	if      (psi.modeMatches(ix,3,mode  )) sign =  1;
        else if (psi.modeMatches(ix,3,modeCC)) sign = -1;
        else  continue;
        const Particle pp     = psi.decayProducts()[ix].at( sign*2212)[0];
        const Particle lamBar = psi.decayProducts()[ix].at(-sign*3122)[0];
        FourMomentum pL = pp.mom()+lamBar.mom();
        double mass = pL.mass();
        _h[0]->fill(mass/GeV);
        double delta = mass-pp.mass()-lamBar.mass();
        _h[1]->fill(delta/GeV);
        if (delta>.15) continue;
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(psi.decaying()[ix].mom().betaVec());
        FourMomentum pProton = boost.transform(pp.mom());
        pL = boost.transform(pL);
        LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pL.betaVec());
        pProton = boost2.transform(pProton);
        const double cosp = -pL.p3().unit().dot(pProton.p3().unit());
        _h[2]->fill(cosp);
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
    const map<PdgId,unsigned int> mode   = { { 2212,1}, {-3122,1}, {-321,1} };
    const map<PdgId,unsigned int> modeCC = { {-2212,1}, { 3122,1}, { 321,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2004_I650570);

}
