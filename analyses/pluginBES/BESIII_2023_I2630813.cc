// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief D_s+ -> omega eta pi+
  class BESIII_2023_I2630813 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2630813);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==431);
      declare(ufs, "UFS");
      DecayedParticles DS(ufs);
      DS.addStable(PID::PI0);
      DS.addStable(PID::K0S);
      DS.addStable(PID::ETA);
      DS.addStable(PID::OMEGA);
      DS.addStable(PID::ETAPRIME);
      declare(DS, "DS");
      // histograms
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles DS = apply<DecayedParticles>(event, "DS");
      // loop over particles
      for (unsigned int ix=0; ix<DS.decaying().size(); ++ix) {
        int sign = 1;
        if (DS.decaying()[ix].pid()>0 && DS.modeMatches(ix,3,mode)) {
          sign=1;
        }
        else if (DS.decaying()[ix].pid()<0 && DS.modeMatches(ix,3,modeCC)) {
          sign=-1;
        }
        else {
          continue;
        }
        const Particle& eta   = DS.decayProducts()[ix].at(     221)[0];
        const Particle& omega = DS.decayProducts()[ix].at(     223)[0];
        const Particle& pip   = DS.decayProducts()[ix].at(sign*211)[0];
        _h[0]->fill((  eta.mom()+pip.mom()).mass());
        _h[1]->fill((omega.mom()+pip.mom()).mass());
        _h[2]->fill((omega.mom()+eta.mom()).mass());
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
    const map<PdgId,unsigned int> mode   = { { 211,1},{ 221,1}, { 223,1}};
    const map<PdgId,unsigned int> modeCC = { {-211,1},{ 221,1}, { 223,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2630813);

}
