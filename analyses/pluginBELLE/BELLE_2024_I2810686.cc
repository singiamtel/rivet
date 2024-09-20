// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B0 -> pi- and B+ rho0 semi-leptonic
  class BELLE_2024_I2810686 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2024_I2810686);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_h[ix],1,1,1+ix);
        book(_c[ix],"TMP/c_"+toString(ix));
      }
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
      // Loop over B0 Mesons
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::B0)) {
        _c[0]->fill();
        if(p.pid()<0) {
          if (isSemileptonicDecay(p, {PID::PIPLUS, PID::ELECTRON, PID::NU_EBAR}) ||
              isSemileptonicDecay(p, {PID::PIPLUS, PID::MUON,     PID::NU_MUBAR})) {
            _h[0]->fill(q2(p, PID::PIPLUS));
          }
        }
        else {
          if (isSemileptonicDecay(p, {PID::PIMINUS, PID::POSITRON, PID::NU_E }) ||
              isSemileptonicDecay(p, {PID::PIMINUS, PID::ANTIMUON, PID::NU_MU})) {
            _h[0]->fill(q2(p, PID::PIMINUS));
          }
        }
      }
      for(const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==PID::BPLUS)) {
        _c[1]->fill();
        if(p.pid()<0) {
          if (isSemileptonicDecay(p, {PID::RHO0, PID::ELECTRON, PID::NU_EBAR}) ||
              isSemileptonicDecay(p, {PID::RHO0, PID::MUON,     PID::NU_MUBAR})) {
            _h[1]->fill(q2(p, PID::RHO0));
          }
        }
        else {
          if (isSemileptonicDecay(p, {PID::RHO0, PID::POSITRON, PID::NU_E }) ||
              isSemileptonicDecay(p, {PID::RHO0, PID::ANTIMUON, PID::NU_MU})) {
            _h[1]->fill(q2(p, PID::RHO0));
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_h[ix], 0.5e4/ *_c[ix]);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    CounterPtr _c[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2024_I2810686);

}
