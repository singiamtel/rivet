// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> J/psi pi+ pi-
  class BESII_2007_I728304 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2007_I728304);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::JPSI);
      declare(psi, "psi");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0; ix<psi.decaying().size(); ++ix) {
        if (!psi.modeMatches(ix,3,mode)) continue;
        const Particle & pip = psi.decayProducts()[ix].at( 211)[0];
        const Particle & pim = psi.decayProducts()[ix].at(-211)[0];
        FourMomentum pSigma = pip.momentum()+pim.momentum();
        _h[0]->fill(pSigma.mass()/GeV);
        const double cSigma = pSigma.z()/pSigma.p3().mod();
        _h[1]->fill(cSigma);
        LorentzTransform boost = LorentzTransform::mkFrameTransformFromBeta(pSigma.betaVec());
        FourMomentum ppip = boost.transform(pip.momentum());
        const double cPi=pSigma.p3().unit().dot(ppip.p3().unit());
        _h[2]->fill(cPi);
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
    const map<PdgId,unsigned int> mode = { { 211,1}, {-211,1}, { 443,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2007_I728304);

}
