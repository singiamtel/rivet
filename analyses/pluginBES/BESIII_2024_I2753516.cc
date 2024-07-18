// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief eta spectra
  class BESIII_2024_I2753516 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2024_I2753516);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(UnstableParticles(Cuts::pid==221), "UFS");
      // histos
      vector<double> energies = { 2.0, 2.2, 2.396, 2.6444, 2.9, 3.05, 3.5000, 3.671};
      unsigned int iloc;
      for(iloc=0;iloc<energies.size();++iloc) {
        if(isCompatibleWithSqrtS(energies[iloc])) break;
      }
      if (iloc==energies.size())
        MSG_ERROR("Beam energy incompatible with analysis.");
      book(_h,1,1,1+iloc);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for(const Particle& p : ufs.particles()) {
        _h->fill(p.momentum().p3().mod());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h,1./sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2024_I2753516);

}
