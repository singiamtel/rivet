// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D+ -> K0S pi+ eta
  class BESIII_2023_I2696635 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2696635);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==411);
      declare(ufs, "UFS");
      DecayedParticles DP(ufs);
      DP.addStable(PID::PI0);
      DP.addStable(PID::K0S);
      DP.addStable(PID::ETA);
      declare(DP,"DP");
      // histograms#
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DP = apply<DecayedParticles>(event, "DP");
      // loop over particles
      for (unsigned int ix=0; ix<DP.decaying().size(); ++ix) {
        int sign = 1;
        if (DP.decaying()[ix].pid()>0 && DP.modeMatches(ix,3,mode)) {
          sign=1;
        }
        else if  (DP.decaying()[ix].pid()<0 && DP.modeMatches(ix,3,modeCC)) {
          sign=-1;
        }
        else {
          continue;
        }
        const Particle& pip = DP.decayProducts()[ix].at( sign*211)[0];
        const Particle& KS0 = DP.decayProducts()[ix].at(      310)[0];
        const Particle& eta = DP.decayProducts()[ix].at(      221)[0];
        _h[0]->fill((pip.mom()+KS0.mom()).mass()/GeV);
        _h[1]->fill((eta.mom()+KS0.mom()).mass()/GeV);
        _h[2]->fill((pip.mom()+eta.mom()).mass()/GeV);
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
    const map<PdgId,unsigned int> mode   = { { 211,1}, { 310,1}, { 221,1}};
    const map<PdgId,unsigned int> modeCC = { {-211,1}, { 310,1}, { 221,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2696635);

}
