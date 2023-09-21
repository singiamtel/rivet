// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  kaon production at low energies
  class BESIII_2021_I1868813 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1868813);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");
      book(_c_kaons   , "/TMP/sigma_kaons");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      unsigned int count = ufs.particles(Cuts::pid==310).size();
      _c_kaons->fill(count);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate1DPtr mult;
      book(mult,1, 1, 1);
      double fact = crossSection()/nanobarn/sumOfWeights();
      double sigma = _c_kaons->val()*fact;
      double error = _c_kaons->err()*fact;
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS()/GeV, b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _c_kaons;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1868813);

}
