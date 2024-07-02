// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > p pbar p nbar pi- +cc
  class BESIII_2022_I2593281 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2022_I2593281);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_sigma, "/TMP/n", refData(1,1,1));
      book(_est, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if ( ntotal == 5 &&
         ((nCount[ 2212]==2 && nCount[-2212]==1 && nCount[-2112]==1 && nCount[-211]==1) ||
          (nCount[-2212]==2 && nCount[ 2212]==1 && nCount[ 2112]==1 && nCount[ 211]==1))) _sigma->fill(sqrtS());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma,crossSection()/ sumOfWeights() /femtobarn);
      barchart(_sigma, _est);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma;
    Estimate1DPtr _est;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2022_I2593281);

}
