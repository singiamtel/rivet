// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief 
  class BESIII_2023_I2670262 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2670262);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable( 3122);
      psi.addStable(-3122);
      psi.addStable( 3112);
      psi.addStable(-3112);
      psi.addStable( 3222);
      psi.addStable(-3222);
      declare(psi, "psi");
      for(unsigned int ix=0;ix<4;++ix) {
        for (unsigned int iy=0; iy<3; ++iy)
          book(_h[ix][iy], 1+ix, 1, 1+iy);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0;ix<psi.decaying().size();++ix) {
        int isgn,ipi,isig=3112,imode;
       	if     (psi.modeMatches(ix,3,mode1)) {
          isgn = -1;
          ipi  =  1;
          imode=  0;
        }
        else if(psi.modeMatches(ix,3,mode2)) {
          isgn =  1;
          ipi  = -1;
          imode=  1;
        }
        else if(psi.modeMatches(ix,3,mode3)) {
          isgn = -1;
          ipi  = -1;
          isig = 3222;
          imode= 2;
        }
        else if(psi.modeMatches(ix,3,mode4)) {
          isgn =  1;
          ipi  =  1;
          isig = 3222;
          imode= 3;
        }
        else continue;
       	const Particle& lam = psi.decayProducts()[ix].at( isgn*3122)[0];
       	const Particle& sig = psi.decayProducts()[ix].at(-isgn*isig)[0];
       	const Particle& pi  = psi.decayProducts()[ix].at( ipi*211  )[0];
        _h[imode][0]->fill((lam.mom()+pi .mom()).mass());
        _h[imode][1]->fill((sig.mom()+pi .mom()).mass());
        _h[imode][2]->fill((lam.mom()+sig.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _h[4][3];
    const map<PdgId,unsigned int> mode1 = { {-3122,1}, { 3112,1}, { 211,1} };
    const map<PdgId,unsigned int> mode2 = { { 3122,1}, {-3112,1}, {-211,1} };
    const map<PdgId,unsigned int> mode3 = { {-3122,1}, { 3222,1}, {-211,1} };
    const map<PdgId,unsigned int> mode4 = { { 3122,1}, {-3222,1}, { 211,1} };
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2670262);

}
