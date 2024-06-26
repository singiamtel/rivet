// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Upsilon(2S) -> Upsilon(1S) pi+pi-
  class ARGUS_1987_I238071 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1987_I238071);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==100553);
      declare(ufs, "UFS");
      DecayedParticles UPS(ufs);
      UPS.addStable( 553);
      declare(UPS, "UPS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_h[ix], 1+ix, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles UPS = apply<DecayedParticles>(event, "UPS");
      // loop over particles
      for (unsigned int ix=0; ix<UPS.decaying().size(); ++ix) {
        // check pi+pi- upslion(1S) decay mode
      	if (!UPS.modeMatches(ix,3,mode)) continue;
        const Particle& ups1 = UPS.decayProducts()[ix].at( 553)[0];
        const Particle& pip  = UPS.decayProducts()[ix].at( 211)[0];
        const Particle& pim  = UPS.decayProducts()[ix].at(-211)[0];
        double mpipi=(pip.momentum()+pim.momentum()).mass();
        _h[0]->fill(mpipi);
        double mpi = pip.mass()+pim.mass();
        double x = (mpipi-mpi)/(UPS.decaying()[ix].mass()-ups1.mass()-mpi);
        _h[1]->fill(x);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    const map<PdgId,unsigned int> mode = { { 553,1},{ 211,1}, {-211,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1987_I238071);

}
