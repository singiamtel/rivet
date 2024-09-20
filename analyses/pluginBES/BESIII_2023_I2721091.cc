// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> pbar Sigma+ K0S
  class BESIII_2023_I2721091 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2721091);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable( 3222);
      psi.addStable(-3222);
      declare(psi, "psi");
      for (unsigned int iy=0; iy<3; ++iy)
        book(_h[iy], 1, 1, 1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0;ix<psi.decaying().size();++ix) {
        int isign=1;
       	if     (psi.modeMatches(ix,3,mode1)) isign =  1;
        else if(psi.modeMatches(ix,3,mode2)) isign = -1;
        else continue;
       	const Particle& k0   = psi.decayProducts()[ix].at( 310       )[0];
       	const Particle& sigp = psi.decayProducts()[ix].at( isign*3222)[0];
       	const Particle& prot = psi.decayProducts()[ix].at(-isign*2212)[0];
        _h[0]->fill((prot.mom()+k0  .mom()).mass());
        _h[1]->fill((sigp.mom()+k0  .mom()).mass());
        _h[2]->fill((sigp.mom()+prot.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode1 = { {-2212,1}, { 3222,1}, { 310,1} };
    const map<PdgId,unsigned int> mode2 = { { 2212,1}, {-3222,1}, { 310,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2721091);

}
