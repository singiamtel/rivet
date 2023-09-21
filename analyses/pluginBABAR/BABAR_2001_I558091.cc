// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief J/psi rate
  class BABAR_2001_I558091 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2001_I558091);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(UnstableParticles(), "UFS");
      book(_c_jpsi, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      double nJsi=0.;
      for (const Particle& p : ufs.particles()) {
        if (p.pid()==443) {
          bool fs = true;
          for (const Particle & child : p.children()) {
            if(child.pid()==443) {
              fs = false;
              break;
            }
          }
          if(fs) nJsi += 1.;
        }
      }
      _c_jpsi->fill(Ecm, nJsi);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      scale(_c_jpsi, fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_jpsi;
    const string Ecm = "10.57";
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(BABAR_2001_I558091);


}
