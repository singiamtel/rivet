// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {

  /// @brief Measurement of forward neutron(+ antineutron) production in proton-proton Collisions at 13 TeV
  class LHCF_2023_I2658888 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LHCF_2023_I2658888);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableFinalState(), "UFS");

      // Book histograms
      book(_h_eta_xf, 1, 1, 1);

    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const UnstableFinalState &ufs = apply<UnstableFinalState> (event, "UFS");

      for (const Particle& p: ufs.particles() ) {
        // select eta meson
        if(p.abspid() != PID::ETA) continue;
        // select pt range
        const double pt = abs(p.pT()/GeV);
        if(pt >= 1.10) continue; // fill considering both LHCf (Arm1 and Arm2) sides
        // compute xf variable
        const double pz = abs(p.pz()/GeV);
        const double ss = sqrtS()/GeV;
        const double xf = 2.*pz/ss;

        _h_eta_xf->fill( xf, xf );
      }

    }

    /// Normalise histograms etc., after the run
    void finalize() {

      const double sumWeights = 2.*sumOfWeights(); // scale considering the LHCf Arm2 side alone

      scale(_h_eta_xf, 1./(sumWeights)); // scale to production rate

    }
    //@}

  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_eta_xf;
    //@}

  };

  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(LHCF_2023_I2658888);

}
