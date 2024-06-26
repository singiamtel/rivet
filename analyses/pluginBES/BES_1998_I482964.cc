// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma K+K-pi0
  class BES_1998_I482964 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1998_I482964);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==443);
      declare(ufs, "UFS");
      DecayedParticles PSI(ufs);
      PSI.addStable(PID::K0S);
      PSI.addStable(PID::PI0);
      declare(PSI, "PSI");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_h[ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      static const map<PdgId,unsigned int> & mode1 = { { 321,1}, { 321,1}, { 111,1}, { 22,1}};
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for (unsigned int ix=0;ix<PSI.decaying().size();++ix) {
        if (!PSI.modeMatches(ix,4,mode1)) continue;
        const Particle& Kp  = PSI.decayProducts()[0].at( 321)[0];
        const Particle& Km  = PSI.decayProducts()[0].at(-321)[0];
        const Particle& pi0 = PSI.decayProducts()[0].at( 111)[0];
        FourMomentum pKK = Kp.mom()+Km.mom();
        _h[0]->fill((pKK+pi0.mom()).mass()/GeV);
        _h[1]->fill((Kp.mom()+pi0.mom()).mass()/GeV);
        _h[1]->fill((Km.mom()+pi0.mom()).mass()/GeV);
        _h[2]->fill(pKK.mass()/GeV);
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


  RIVET_DECLARE_PLUGIN(BES_1998_I482964);

}
