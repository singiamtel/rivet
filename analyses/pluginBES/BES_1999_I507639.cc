// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief J/psi -> gamma 2pi+ 2pi-
  class BES_1999_I507639 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BES_1999_I507639);


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
      book(_h,1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the J/psi decays
      DecayedParticles PSI = apply<DecayedParticles>(event, "PSI");
      // loop over particles
      for (unsigned int ix=0; ix<PSI.decaying().size(); ++ix) {
        if (!PSI.modeMatches(ix,5,mode)) continue;
        const Particles& pip  = PSI.decayProducts()[0].at( 211);
        const Particles& pim  = PSI.decayProducts()[0].at(-211);
        FourMomentum ptotal;
        for (unsigned int iy=0; iy<2; ++iy) {
          ptotal+=pip[ix].momentum()+pim[ix].momentum();
        }
        _h->fill(ptotal.mass());
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
    const map<PdgId,unsigned int> mode = { { 211,2}, {-211,2}, { 22,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BES_1999_I507639);

}
