// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief tau -> 3 pi
  class ARGUS_1992_I339268 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1992_I339268);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(111);
      declare(TAU, "TAU");
      // histos
      book(_h_mass,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
      	int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
      	if (!(TAU.modeMatches(ix,4,mode) || TAU.modeMatches(ix,4,modeCC))) {
          continue;
        }
      	const Particles& pip = TAU.decayProducts()[ix].at( 211*sign);
        const Particles & pim = TAU.decayProducts()[ix].at(-211*sign);
        double Q2 = (pip[0].mom()+pim[0].mom()+pim[1].mom()).mass2();
        _h_mass->fill(sqrt(Q2));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h_mass,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass;
    const map<PdgId,unsigned int> mode   = { {-211,2}, { 211,1}, { 16,1}};
    const map<PdgId,unsigned int> modeCC = { { 211,2}, {-211,1}, {-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1992_I339268);

}
