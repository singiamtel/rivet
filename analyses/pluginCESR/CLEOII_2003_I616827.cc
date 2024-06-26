// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 -> pi/rho ell+ nu_ell
  class CLEOII_2003_I616827 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_2003_I616827);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::abspid==511), "UFS");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
      book(_c,"TMP/nB");
    }

    // Calculate the Q2 using mother and daugher meson
    double q2(const Particle& B, int mesonID) {
      FourMomentum q = B.mom() - select(B.children(), Cuts::pid==mesonID)[0];
      return q*q;
    }

    // Check for explicit decay into pdgids
    bool isSemileptonicDecay(const Particle& mother, vector<int> ids) {
      // Trivial check to ignore any other decays but the one in question modulo photons
      const Particles children = mother.children(Cuts::pid!=PID::PHOTON);
      if (children.size()!=ids.size()) return false;
      // Check for the explicit decay
      return all(ids, [&](int i){return count(children, hasPID(i))==1;});
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Loop over B0bar Mesons
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::pid==PID::B0BAR)) {
        _c->fill();
        if (isSemileptonicDecay(p, {PID::PIPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::PIPLUS, PID::MUON,     PID::NU_MUBAR})) {
          _h[0]->fill(q2(p, PID::PIPLUS));
        }
        if (isSemileptonicDecay(p, {PID::RHOPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
            isSemileptonicDecay(p, {PID::RHOPLUS, PID::MUON,     PID::NU_MUBAR})) {
          _h[1]->fill(q2(p, PID::RHOPLUS));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // 1e4 for q0-4 in BR and 0.5 as both e/mu modes
      scale(_h, 0.5*1e4/ *_c);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_2003_I616827);

}
