// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> K- Lambda9 Xibar+
  class BESIII_2024_I2691955 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2691955);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable( 3122);
      psi.addStable(-3122);
      psi.addStable( 3312);
      psi.addStable(-3312);
      declare(psi, "psi");
      for (unsigned int iy=0; iy<3; ++iy)
        book(_h[iy], 1, 1, 1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0;ix<psi.decaying().size();++ix) {
        int isgn=1;
       	if     (psi.modeMatches(ix,3,mode1)) isgn= 1;
        else if(psi.modeMatches(ix,3,mode2)) isgn=-1;
        else continue;
       	const Particle& Km  = psi.decayProducts()[ix].at(-isgn* 321)[0];
       	const Particle& lam = psi.decayProducts()[ix].at( isgn*3122)[0];
       	const Particle& xi  = psi.decayProducts()[ix].at(-isgn*3312)[0];
        _h[0]->fill((Km .mom()+xi .mom()).mass());
        _h[1]->fill((lam.mom()+xi .mom()).mass());
        _h[2]->fill((Km .mom()+lam.mom()).mass());
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
    const map<PdgId,unsigned int> mode1 = { {-321,1}, { 3122,1}, {-3312,1} };
    const map<PdgId,unsigned int> mode2 = { { 321,1}, {-3122,1}, { 3312,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2691955);

}
