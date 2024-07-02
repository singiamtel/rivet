// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  tau -> 6 pion final states
  class CLEO_2001_I535016 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_2001_I535016);


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
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0; ix<TAU.decaying().size(); ++ix) {
        int sign = TAU.decaying()[ix].pid()>0 ? 1 : -1;
        int imode=-1;
        if (TAU.modeMatches(ix,7,mode1  ) ||
            TAU.modeMatches(ix,7,mode1CC)) imode = 0;
        else if (TAU.modeMatches(ix,7,mode2) ||
                 TAU.modeMatches(ix,7,mode2CC)) imode = 1;
        else continue;
        const Particles& pip = TAU.decayProducts()[ix].at( 211*sign);
        const Particles& pim = TAU.decayProducts()[ix].at(-211*sign);
        const Particles& pi0 = TAU.decayProducts()[ix].at( 111);
        FourMomentum ptotal;
        for (unsigned int ix=0;ix<pim.size();++ix) ptotal +=pim[ix].momentum();
        for (unsigned int ix=0;ix<pip.size();++ix) ptotal +=pip[ix].momentum();
        for (unsigned int ix=0;ix<pi0.size();++ix) ptotal +=pi0[ix].momentum();
        _h[imode]->fill(ptotal.mass());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode1   = { {-211,2}, { 211,1},{111,3},{ 16,1}};
    const map<PdgId,unsigned int> mode1CC = { { 211,2}, {-211,1},{111,3},{-16,1}};
    const map<PdgId,unsigned int> mode2   = { {-211,3}, { 211,2},{111,1},{ 16,1}};
    const map<PdgId,unsigned int> mode2CC = { { 211,3}, {-211,2},{111,1},{-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_2001_I535016);

}
