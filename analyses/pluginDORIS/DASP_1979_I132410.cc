// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief inclusive eta production
  class DASP_1979_I132410 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DASP_1979_I132410);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_eta, "/TMP/neta");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      _c_eta->fill(ufs.particles(Cuts::pid==221).size());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c_eta,1./sumOfWeights()*crossSection()/nanobarn);
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(_c_eta->val(), _c_eta->err());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_eta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DASP_1979_I132410);


}
