// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> KS0 pi+ pi0 pi0
  class BESIII_2023_I2662580 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2662580);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DD(ufs);
      DD.addStable(PID::PI0);
      DD.addStable(PID::K0S);
      declare(DD, "DD");
      // histograms
      for (unsigned int ix=0; ix<7; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DD = apply<DecayedParticles>(event, "DD");
      // loop over particles
      for (unsigned int ix=0; ix<DD.decaying().size(); ++ix) {
        int sign = DD.decaying()[ix].pid()/DD.decaying()[ix].abspid();
        if (!DD.modeMatches(ix,4,mode) && !DD.modeMatches(ix,4,modeCC)) continue;
        const Particle & KS0 = DD.decayProducts()[ix].at(      310)[0];
        const Particle & pip = DD.decayProducts()[ix].at( sign*211)[0];
        const Particles& pi0 = DD.decayProducts()[ix].at(      111);
        const double mpipi = (pi0[0].mom()+pi0[1].mom()).mass();
        if (mpipi>0.46 && mpipi<0.52) continue;
        _h[3]->fill(mpipi);
        _h[2]->fill((KS0.mom()+pip.mom()).mass());
        _h[5]->fill((KS0.mom()+pi0[0].mom()+pi0[1].mom()).mass());
        _h[6]->fill((pip.mom()+pi0[0].mom()+pi0[1].mom()).mass());
        for(unsigned int ix=0; ix<2; ++ix) {
          _h[0]->fill((KS0.mom()+pi0[ix].mom()).mass());
          _h[1]->fill((pip.mom()+pi0[ix].mom()).mass());
          _h[4]->fill((KS0.mom()+pip.mom()+pi0[ix].mom()).mass());
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
	    normalize(_h);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[7];
    const map<PdgId,unsigned int> mode   = { {310,1}, { 211,1}, { 111,2} };
    const map<PdgId,unsigned int> modeCC = { {310,1}, {-211,1}, { 111,2} };
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2662580);

}
