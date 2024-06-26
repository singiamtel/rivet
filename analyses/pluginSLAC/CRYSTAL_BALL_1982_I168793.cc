// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Diphoton to hadrons at centre-of-mass energies between 0.72 and 1.8 GeV
  class CRYSTAL_BALL_1982_I168793 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1982_I168793);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if (inRange(sqrtS()/GeV,0.72,1.8)) {
        book(_mult, "TMP/mult", refData(1, 1, 1));
        if (inRange(sqrtS()/GeV,1.04,1.48)) book(_h_cTheta,2,1,1);
      }
      else
       	throw Error("Invalid CMS energy for CRYSTAL_BALL_1990_I294492");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      for (const Particle & p : part) {
        if (p.pid()!=PID::PI0) vetoEvent;
      }
      const double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if (cTheta<=0.7) _mult->fill(sqrtS()/GeV);
      if (_h_cTheta ) _h_cTheta->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta) scale(_h_cTheta,fact);
      scale(_mult, fact);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_mult,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _mult, _h_cTheta;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1982_I168793);

}
