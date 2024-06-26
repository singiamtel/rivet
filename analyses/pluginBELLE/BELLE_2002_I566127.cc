// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Lambda_c decays
  class BELLE_2002_I566127 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2002_I566127);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==4122);
      declare(ufs, "UFS");
      DecayedParticles LAMBDAC(ufs);
      LAMBDAC.addStable(PID::PI0);
      LAMBDAC.addStable(PID::K0S);
      LAMBDAC.addStable( 3122);
      LAMBDAC.addStable(-3122);
      LAMBDAC.addStable( 3222);
      LAMBDAC.addStable(-3222);
      declare(LAMBDAC, "LAMBDAC");
      // histos
      for (unsigned int ix=0;ix<3;++ix) {
        book(_h[ix],1+ix,1,1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles LAMBDAC = apply<DecayedParticles>(event, "LAMBDAC");
      // loop over particles
      for (unsigned int ix=0;ix<LAMBDAC.decaying().size();++ix) {
        int sign = LAMBDAC.decaying()[ix].pid()/LAMBDAC.decaying()[ix].abspid();
        if ((LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode1  )) ||
            (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,mode1CC))) {
          const Particle& sig = LAMBDAC.decayProducts()[ix].at( sign*3222)[0];
          const Particle& Kp  = LAMBDAC.decayProducts()[ix].at( sign*321 )[0];
          const Particle& Km  = LAMBDAC.decayProducts()[ix].at(-sign*321 )[0];
          _h[0]->fill((Kp .mom()+Km.mom()).mass());
          _h[1]->fill((sig.mom()+Km.mom()).mass());
        }
        else if ((LAMBDAC.decaying()[ix].pid()>0 && LAMBDAC.modeMatches(ix,3,mode2)) ||
                 (LAMBDAC.decaying()[ix].pid()<0 && LAMBDAC.modeMatches(ix,3,mode2CC))) {
          const Particle& lam = LAMBDAC.decayProducts()[ix].at(sign*3122)[0];
          const Particle& K0  = LAMBDAC.decayProducts()[ix].at(     310 )[0];
          _h[2]->fill((lam.mom()+K0.mom()).mass());
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
    const map<PdgId,unsigned int> mode1   = { { 3222,1}, { 321,1}, {-321,1} };
    const map<PdgId,unsigned int> mode1CC = { {-3222,1}, { 321,1}, {-321,1} };
    const map<PdgId,unsigned int> mode2   = { { 3122,1}, { 310,1}, { 321,1} };
    const map<PdgId,unsigned int> mode2CC = { {-3122,1}, { 310,1}, {-321,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2002_I566127);

}
