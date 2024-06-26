// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> K- pi+ pi- nu_tau
  class CLEO_2000_I525698 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2000_I525698);


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
      declare(TAU, "TAU");
      // histogram
      for (unsigned int ix=0;ix<3;++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0;ix<TAU.decaying().size();++ix) {
        int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
        if (!(TAU.modeMatches(ix,4,mode  ) || TAU.modeMatches(ix,4,modeCC))) continue;
        const Particle & pip = TAU.decayProducts()[ix].at( 211*sign)[0];
        const Particle & pim = TAU.decayProducts()[ix].at(-211*sign)[0];
        const Particle & Km  = TAU.decayProducts()[ix].at(-321*sign)[0];
        _h[0]->fill((Km.mom()+pip.mom()+pim.mom()).mass());
        _h[1]->fill((Km.mom()+pip.mom()).mass());
        _h[2]->fill((pip.mom()+pim.mom()).mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
	     normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode   = { {-321,1}, { 211,1},{-211,1},{ 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 321,1}, {-211,1},{ 211,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2000_I525698);

}
