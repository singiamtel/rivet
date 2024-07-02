// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K+K-pi+pi0
  class BESIII_2020_I1808166 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1808166);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");
      // histograms
      for (unsigned int ix=0; ix<10; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for (unsigned int ix=0; ix<DD.decaying().size(); ++ix) {
      	int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
      	if ( !DD.modeMatches(ix,4,mode) && !DD.modeMatches(ix,4,modeCC)) continue;
          const Particle & Kp  = DD.decayProducts()[ix].at( sign*321)[0];
          const Particle & Km  = DD.decayProducts()[ix].at(-sign*321)[0];
          const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
          const Particle & pi0 = DD.decayProducts()[ix].at(      111)[0];
          _h[0]->fill((Kp .mom()+Km .mom()).mass());
          _h[1]->fill((Kp .mom()+pip.mom()).mass());
          _h[2]->fill((Kp .mom()+pi0.mom()).mass());
          _h[3]->fill((Km .mom()+pip.mom()).mass());
          _h[4]->fill((Km .mom()+pi0.mom()).mass());
          _h[5]->fill((pip.mom()+pi0.mom()).mass());
          _h[6]->fill((Kp .mom()+Km .mom()+pip.mom()).mass());
          _h[7]->fill((Kp .mom()+Km .mom()+pi0.mom()).mass());
          _h[8]->fill((Kp .mom()+pip.mom()+pi0.mom()).mass());
          _h[9]->fill((Km .mom()+pip.mom()+pi0.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[10];
    const map<PdgId,unsigned int> mode   = { {321,1}, {-321,1}, { 211,1}, { 111,1} };
    const map<PdgId,unsigned int> modeCC = { {321,1}, {-321,1}, {-211,1}, { 111,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1808166);

}
