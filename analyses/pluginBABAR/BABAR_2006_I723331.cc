// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B0 -> D*+ D*- KS0
  class BABAR_2006_I723331 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I723331);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::PI0);
      BB.addStable( 413);
      BB.addStable(-413);
      BB.addStable( PID::K0S);
      declare(BB, "BB");
      // histograms
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode = { { 413,1}, { -413,1}, { 310,1}};
      // loop over particles
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	if ( !BB.modeMatches(ix,3,mode) ) continue;
	const Particle & K0  = BB.decayProducts()[ix].at( 310)[0];
	const Particle & DsM = BB.decayProducts()[ix].at( 413)[0];
	const Particle & DsP = BB.decayProducts()[ix].at(-413)[0];
	_h->fill((K0.momentum()+DsM.momentum()).mass());
	_h->fill((K0.momentum()+DsP.momentum()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I723331);

}
