// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  chi_c1,2 -> Lambda Lambdabar0 phi
  class BESIII_2024_I2792993 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2792993);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::pid==20443 or
						Cuts::pid==445);
      declare(ufs, "UFS");
      DecayedParticles chi(ufs);
      chi.addStable( PID::PHI);
      chi.addStable( PID::LAMBDA);
      chi.addStable(-PID::LAMBDA);
      declare(chi, "chi");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
        for(unsigned int iy=0;iy<3;++iy)
          book(_h[ix][iy],1+ix,1,1+iy);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 3122,1}, {-3122,1}, { 333,1} };
      DecayedParticles chi = apply<DecayedParticles>(event, "chi");
      for(unsigned int ix=0;ix<chi.decaying().size();++ix) {
	if(!chi.modeMatches(ix,3,mode)) continue;
        unsigned int imode = chi.decaying()[ix].pid()==20443 ? 0 : 1;
	const Particle & lam = chi.decayProducts()[ix].at( 3122)[0];
	const Particle & lamb = chi.decayProducts()[ix].at(-3122)[0];
	const Particle & phi  = chi.decayProducts()[ix].at(  333)[0];
	_h[imode][2]->fill((lam .momentum()+lamb.momentum()).mass());
	_h[imode][0]->fill((lam .momentum()+phi .momentum()).mass());
	_h[imode][1]->fill((lamb.momentum()+phi .momentum()).mass());
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
        for(unsigned int iy=0;iy<3;++iy)
          normalize(_h[ix][iy],1.,false);
    }

    /// @}

    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][3];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2792993);

}
