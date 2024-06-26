// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> J/psi eta' K+
  class LHCB_2023_I2643022 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCB_2023_I2643022);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::ETAPRIME);
      BP.addStable(PID::JPSI);
      declare(BP, "BP");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> mode   = { { 321,1}, { 331,1},{ 443,1}};
      static const map<PdgId,unsigned int> modeCC = { {-321,1}, { 331,1},{ 443,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for (unsigned int ix=0; ix<BP.decaying().size(); ++ix) {
      	int sign = 1;
      	if (BP.decaying()[ix].pid()>0 && BP.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if (BP.decaying()[ix].pid()<0 && BP.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else {
      	  continue;
        }
        const Particle& Kp  = BP.decayProducts()[ix].at( sign*321)[0];
        const Particle& eta = BP.decayProducts()[ix].at(      331)[0];
        const Particle& psi = BP.decayProducts()[ix].at(      443)[0];
        _h[0]->fill((Kp .mom()+eta.mom()).mass());
        _h[1]->fill((psi.mom()+eta.mom()).mass());
        _h[2]->fill((psi.mom()+Kp .mom()).mass());
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
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LHCB_2023_I2643022);

}
