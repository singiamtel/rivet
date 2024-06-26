// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief chi_c2 decays
  class BESIII_2017_I1505590 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1505590);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PI0);
      chi.addStable( PID::K0S);
      declare(chi, "chi");
      // histograms
      for(unsigned int ix=0; ix<4; ++ix) {
        if (ix<2)  book(_h[ix], 1, 1, 1+ix);
        book(_h[2+ix], 2, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      // loop over particles
      for (unsigned int ix=0; ix<chi.decaying().size(); ++ix) {
        int sign=1;
        // K+ K- pi0
        if (chi.modeMatches(ix,3,mode1)) {
          const Particle& Kp  = chi.decayProducts()[ix].at( 321)[0];
          const Particle& Km  = chi.decayProducts()[ix].at(-321)[0];
          const Particle& pi0 = chi.decayProducts()[ix].at( 111)[0];
          _h[0]->fill((Kp.mom()+pi0.mom()).mass()/GeV);
          _h[0]->fill((Km.mom()+pi0.mom()).mass()/GeV);
          _h[1]->fill((Kp.mom()+Km .mom()).mass()/GeV);
          continue;
        }
        // pi+ pi- pi0
        else if (chi.modeMatches(ix,3,mode3)) {
          const Particle & pip  = chi.decayProducts()[ix].at( 211)[0];
          const Particle & pim  = chi.decayProducts()[ix].at(-211)[0];
          const Particle & pi0 = chi.decayProducts()[ix].at( 111)[0];
          _h[5]->fill((pip.mom()+pi0.mom()).mass()/GeV);
          _h[5]->fill((pim.mom()+pi0.mom()).mass()/GeV);
          continue;
        }
        else if (chi.modeMatches(ix,3,mode2)) {
          sign =  1;
        }
        else if (chi.modeMatches(ix,3,mode2CC)) {
          sign = -1;
        }
        else {
          continue;
        }
        const Particle& KS0 = chi.decayProducts()[ix].at(      310)[0];
        const Particle& Kp  = chi.decayProducts()[ix].at( sign*321)[0];
        const Particle & pim = chi.decayProducts()[ix].at(-sign*211)[0];
        _h[2]->fill((Kp .mom()+pim.mom()).mass()/GeV);
        _h[3]->fill((KS0.mom()+pim.mom()).mass()/GeV);
        _h[4]->fill((KS0.mom()+Kp .mom()).mass()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[6];
    const map<PdgId,unsigned int> mode1   = { { 321,1}, {-321,1}, { 111,1} };
    const map<PdgId,unsigned int> mode2   = { { 310,1}, { 321,1}, {-211,1} };
    const map<PdgId,unsigned int> mode2CC = { { 310,1}, {-321,1}, { 211,1} };
    const map<PdgId,unsigned int> mode3   = { { 211,1}, {-211,1}, { 111,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1505590);

}
