// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief psi(2S) -> Lambda Lambdabar eta'
  class BESIII_2023_I2702517 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2702517);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==100443);
      declare(ufs, "UFS");
      DecayedParticles psi(ufs);
      psi.addStable(PID::PI0);
      psi.addStable(PID::K0S);
      psi.addStable(PID::ETA);
      psi.addStable(PID::ETAPRIME);
      psi.addStable(PID::OMEGA);
      psi.addStable( 3122);
      psi.addStable(-3122);
      declare(psi, "psi");
      for (unsigned int ix=0; ix<2; ++ix) {
        for (unsigned int iy=0; iy<2; ++iy) {
          book(_h[ix][iy], 1+ix, 1, 1+iy);
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles psi = apply<DecayedParticles>(event, "psi");
      // loop over particles
      for (unsigned int ix=0;ix<psi.decaying().size();++ix) {
       	if (!psi.modeMatches(ix,3,mode)) continue;
       	const Particle& eta  = psi.decayProducts()[ix].at( 331)[0];
       	const Particle& lam  = psi.decayProducts()[ix].at( 3122)[0];
       	const Particle& lbar = psi.decayProducts()[ix].at(-3122)[0];
        for (unsigned int iy=0; iy<2; ++iy) {
          _h[0][iy]->fill((lbar.mom()+eta .mom()).mass());
          _h[0][iy]->fill((lam .mom()+eta .mom()).mass());
          _h[1][iy]->fill((lam .mom()+lbar.mom()).mass());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<2;++ix) {
        normalize(_h[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    const map<PdgId,unsigned int> mode = { { 3122,1}, {-3122,1}, { 331,1} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2702517);

}
