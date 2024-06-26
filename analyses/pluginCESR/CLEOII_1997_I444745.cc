// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> pi+ pi- pi- eta nu_tau
  class CLEOII_1997_I444745 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOII_1997_I444745);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      TAU.addStable(221);
      declare(TAU, "TAU");
      // histogram
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for(unsigned int ix=0;ix<TAU.decaying().size();++ix) {
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
      	if ( !(TAU.modeMatches(ix,5,mode) || TAU.modeMatches(ix,5,modeCC)) ) continue;
      	const Particle&  pip = TAU.decayProducts()[ix].at( 211*sign)[0];
      	const Particles& pim = TAU.decayProducts()[ix].at(-211*sign);
        const Particle&  eta = TAU.decayProducts()[ix].at( 221)[0];
        FourMomentum phad = pip.mom()+eta.mom();
        bool hasF1=false;
        for (unsigned int iy=0; iy<2; ++iy) {
          double mf1 = (phad+pim[iy].mom()).mass();
          hasF1 |= abs(mf1-1.2819) < 0.036;
        }
        if (hasF1) continue;
        phad += pim[0].mom() + pim[1].mom();
        _h->fill(phad.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    const map<PdgId,unsigned int> mode   = { { 221,1}, { 211,1},{-211,2},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 221,1}, {-211,1},{ 211,2},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOII_1997_I444745);

}
