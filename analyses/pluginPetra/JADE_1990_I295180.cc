// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi0pi0
  class JADE_1990_I295180 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1990_I295180);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if (inRange(sqrtS()/GeV,2.,3.5)) {
        book(_npipi, "TMP/pipi", refData(1, 1, 1));
      }
      else {
       	throw Error("Invalid CMS energy for CRYSTAL_BALL_1990_I294492");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      for (const Particle& p : part) {
        if (p.pid()!=PID::PI0) vetoEvent;
      }
      const double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if (cTheta<=0.3) _npipi->fill(sqrtS()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npipi, crossSection()/nanobarn/sumOfWeights());
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_npipi,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _npipi;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(JADE_1990_I295180);

}
