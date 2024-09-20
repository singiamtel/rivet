// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Sigma+ Sigmabar- + omega/phi
  class BESIII_2023_I2713388 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2713388);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::PHI);
      psi.addStable(PID::OMEGA);
      psi.addStable( 3222);
      psi.addStable(-3222);
      declare(psi, "psi");
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<3; ++iy) {
          book(_h[ix][iy], 1+ix, 1, 1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0;ix<psi.decaying().size();++ix) {
        unsigned int iloc=0;
       	if     (psi.modeMatches(ix,3,mode1)) iloc=0;
        else if(psi.modeMatches(ix,3,mode2)) iloc=1;
        else continue;
       	const Particle& vecm = psi.decayProducts()[ix].at(iloc==0 ? 223 : 333)[0];
       	const Particle& sigp = psi.decayProducts()[ix].at( 3222)[0];
       	const Particle& sigm = psi.decayProducts()[ix].at(-3222)[0];
        _h[iloc][0]->fill((sigp.mom()+vecm.mom()).mass());
        _h[iloc][1]->fill((sigm.mom()+vecm.mom()).mass());
        _h[iloc][2]->fill((sigp.mom()+sigm.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<2;++ix) normalize(_h[ix]);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    const map<PdgId,unsigned int> mode1 = { { 3222,1}, {-3222,1}, { 223,1} };
    const map<PdgId,unsigned int> mode2 = { { 3222,1}, {-3222,1}, { 333,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2713388);

}
