// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B- -> Lambda_c+ Lambdabar_c- K-
  class BABAR_2008_I766290 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I766290);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==PID::BPLUS);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 4122);
      BP.addStable(-4122);
      declare(BP,"BP");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 4122,1}, {-4122,1}, { 321,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 4122,1}, {-4122,1}, {-321,1} };
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	// select right decay mode
	int sign = 1;
	if      (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode  )) sign = 1;
	else if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) sign =-1;
	else continue;
	const Particle & lam  = BP.decayProducts()[ix].at( 4122*sign)[0];
	const Particle & lbar = BP.decayProducts()[ix].at(-4122*sign)[0];
	const Particle & Kp   = BP.decayProducts()[ix].at(  321*sign)[0];
	_h[0]->fill((lbar.momentum()+Kp  .momentum()).mass());
	_h[1]->fill((lam .momentum()+lbar.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BABAR_2008_I766290);

}
