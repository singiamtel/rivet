// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief OPAL b-fragmentation measurement for weak B-hadron decays
  /// @author Simone Amoroso
  class OPAL_2003_I599181 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(OPAL_2003_I599181);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_estXbweak, 1, 1, 1);
      book(_histXbweak, "/TMP/Xbweak", refData(1, 1, 1));

      book(_histMeanXbweak, 2, 1, 1);

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      // Get Bottom hadrons
      const Particles bhads = select(ufs.particles(), isBottomHadron);

      for (const Particle& bhad : bhads) {
        // Check for weak decay, i.e. no more bottom present in children
        if (bhad.isLastWith(hasBottom)) {
          const double xp = bhad.E()/meanBeamMom;
          _histXbweak->fill(xp);
          _histMeanXbweak->fill(_histMeanXbweak->bin(1).xMid(), xp);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_histXbweak);
      barchart(_histXbweak, _estXbweak);
    }

    /// @}


  private:

    Histo1DPtr _histXbweak;
    Estimate1DPtr _estXbweak;
    Profile1DPtr _histMeanXbweak;

  };



  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(OPAL_2003_I599181);


}
