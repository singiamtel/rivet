// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"


namespace Rivet {


  /// @brief tau -> hadrons with KS0
  class CLEOII_1996_I415408 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1996_I415408);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      declare(TAU, "TAU");
      // histos
      for (unsigned int ix=0; ix<4; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
        if(ix > 1) continue;
        book(_h[4+ix], 5, 1, 1+ix);
        book(_h[6+ix], 6+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
        int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
        if ((TAU.modeMatches(ix,3,mode1  ) || TAU.modeMatches(ix,3,mode1CC))) {
          const Particle & pim = TAU.decayProducts()[ix].at(-211*sign)[0];
          const Particle & K0  = TAU.decayProducts()[ix].at( 310     )[0];
          _h[0]->fill((pim.mom()+K0.mom()).mass());
        }
        else if (TAU.modeMatches(ix,3,mode2  ) || TAU.modeMatches(ix,3,mode2CC)) {
          const Particle & Km = TAU.decayProducts()[ix].at(-321*sign)[0];
          const Particle & K0 = TAU.decayProducts()[ix].at( 310     )[0];
          _h[1]->fill((Km.mom()+K0.mom()).mass());
        }
        else if (TAU.modeMatches(ix,4,mode3  ) || TAU.modeMatches(ix,4,mode3CC)) {
          const Particle & pim = TAU.decayProducts()[ix].at(-211*sign)[0];
          const Particle & K0  = TAU.decayProducts()[ix].at( 310     )[0];
          const Particle & pi0 = TAU.decayProducts()[ix].at( 111     )[0];
          _h[2]->fill((pim.mom()+K0.mom()+pi0.mom()).mass());
        }
        else if (TAU.modeMatches(ix,4,mode4  ) || TAU.modeMatches(ix,4,mode4CC)) {
          const Particle & Km = TAU.decayProducts()[ix].at(-321*sign)[0];
          const Particle & K0  = TAU.decayProducts()[ix].at( 310     )[0];
          const Particle & pi0 = TAU.decayProducts()[ix].at( 111     )[0];
          _h[3]->fill((Km.mom()+K0.mom()+pi0.mom()).mass());
          _h[4]->fill((Km.mom()+pi0.mom()).mass());
          _h[5]->fill((K0.mom()+pi0.mom()).mass());
        }
        else if (TAU.modeMatches(ix,4,mode5  ) || TAU.modeMatches(ix,4,mode5CC)) {
          const Particle  & pim = TAU.decayProducts()[ix].at(-211*sign)[0];
          const Particles & K0  = TAU.decayProducts()[ix].at( 310     );
          _h[6]->fill((pim.mom()+K0[0].mom()+K0[1].mom()).mass());
          _h[7]->fill((pim.mom()+K0[0].mom()).mass());
          _h[7]->fill((pim.mom()+K0[1].mom()).mass());
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
    Histo1DPtr _h[8];
    const map<PdgId,unsigned int> mode1   = { { 310,1},{-211,1},{ 16,1}};
    const map<PdgId,unsigned int> mode1CC = { { 310,1},{ 211,1},{-16,1}};
    const map<PdgId,unsigned int> mode2   = { { 310,1},{-321,1},{ 16,1}};
    const map<PdgId,unsigned int> mode2CC = { { 310,1},{ 321,1},{-16,1}};
    const map<PdgId,unsigned int> mode3   = { { 310,1},{-211,1},{ 111,1},{ 16,1}};
    const map<PdgId,unsigned int> mode3CC = { { 310,1},{ 211,1},{ 111,1},{-16,1}};
    const map<PdgId,unsigned int> mode4   = { { 310,1},{-321,1},{ 111,1},{ 16,1}};
    const map<PdgId,unsigned int> mode4CC = { { 310,1},{ 321,1},{ 111,1},{-16,1}};
    const map<PdgId,unsigned int> mode5   = { { 310,2},{-211,1},{ 16,1}};
    const map<PdgId,unsigned int> mode5CC = { { 310,2},{ 211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1996_I415408);

}
