// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> J/psi omega K
  class BABAR_2008_I767506 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I767506);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(310);
      BB.addStable(223);
      BB.addStable(443);
      declare(BB, "BB");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 443,1},{ 223,1}, { 321,1}};
      static const map<PdgId,unsigned int> & mode1CC = { { 443,1},{ 223,1}, {-321,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 443,1},{ 223,1}, { 130,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 443,1},{ 223,1}, { 310,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
	unsigned int imode=0;
      	if (BB.modeMatches(ix,3,mode1) ||
	    BB.modeMatches(ix,3,mode1CC))
	  imode=0;
	else if(BB.modeMatches(ix,3,mode2) ||
		BB.modeMatches(ix,3,mode2CC))
	  imode=1;
	else
	  continue;
	const Particle & omega = BB.decayProducts()[ix].at(223)[0];
	const Particle & psi   = BB.decayProducts()[ix].at(443)[0];
	_h[imode]->fill((omega.momentum()+psi.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I767506);

}
