// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class JADE_1986_I231554 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1986_I231554);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,2.,2.6)) {
       	book(_cProton, "TMP/proton", refData(1,1,1));
        book(_h_cTheta,2,1,1);
      }
      else {
        throw Error("Invalid CMS energy for JADE_1986_I231554");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PROTON) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::ANTIPROTON) foundM=true;
      }
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=0.6)    _cProton->fill(sqrtS()/GeV);
      if (_h_cTheta) _h_cTheta->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta) scale(_h_cTheta, fact);
      scale(_cProton, fact);
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_cProton,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cTheta, _cProton;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(JADE_1986_I231554);

}
