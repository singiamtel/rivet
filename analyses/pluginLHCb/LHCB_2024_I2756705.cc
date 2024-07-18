// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B_c -> J/psi pi+pi0
  class LHCB_2024_I2756705 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2024_I2756705);

    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==541);
      declare(ufs, "UFS");
      DecayedParticles BC(ufs);
      BC.addStable( PID::PI0);
      BC.addStable( PID::JPSI);
      declare(BC, "BC");
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 443,1}, { 211,1}, { 111,1} };
      static const map<PdgId,unsigned int> & modeCC = { { 443,1}, {-211,1}, { 111,1} };
      DecayedParticles BC = apply<DecayedParticles>(event, "BC");
      // loop over particles
      for(unsigned int ix=0;ix<BC.decaying().size();++ix) {
	int sign = BC.decaying()[ix].pid()/BC.decaying()[ix].abspid();
	if ((sign== 1 && BC.modeMatches(ix,3,mode  )) ||
	    (sign==-1 && BC.modeMatches(ix,3,modeCC))) {
	  const Particle & pi0 = BC.decayProducts()[ix].at(      111)[0];
	  const Particle & pip = BC.decayProducts()[ix].at( sign*211)[0];
	  FourMomentum ptotal=pi0.momentum()+pip.momentum();
          _h->fill(ptotal.mass());
        }
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


  RIVET_DECLARE_PLUGIN(LHCB_2024_I2756705);

}
