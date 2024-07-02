// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> eta K gamma
  class BABAR_2009_I785439 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2009_I785439);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::ETA);
      BB.addStable(PID::K0S);
      declare(BB, "BB");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 221,1},{ 310,1}, {22,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 221,1},{ 321,1}, {22,1}};
      static const map<PdgId,unsigned int> & mode2CC = { { 221,1},{-321,1}, {22,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
       	int imode=0,sign = 1;
	if(BB.decaying()[ix].abspid()==511) {
	  imode=1;
	  if (!BB.modeMatches(ix,3,mode1)) continue;
	}
	else {
	  imode=0;
	  if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	    sign=1;
	  }
	  else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
	    sign=-1;
	  }
	  else
	    continue;
	}
	const Particle & KK    = BB.decayProducts()[ix].at(imode==0 ?  sign*321 : 310)[0];
      	const Particle & eta   = BB.decayProducts()[ix].at(      221)[0];
	_h[imode]->fill((KK.momentum()+eta.momentum()).mass());
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


  RIVET_DECLARE_PLUGIN(BABAR_2009_I785439);

}
