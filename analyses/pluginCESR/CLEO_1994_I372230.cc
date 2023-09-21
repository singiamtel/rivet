// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-/K+ K-
  class CLEO_1994_I372230 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1994_I372230);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (sqrtS()<1.5*GeV || sqrtS()>5*GeV)
        throw Error("Invalid CMS energy for CLEO_1994_I372230");
      book(_cPiK, "/TMP/nPiK_");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      if (part[0].pid()!=-part[1].pid()) vetoEvent;
      double cTheta(0.);
      bool foundPi(false),foundK(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PIPLUS) {
          foundPi=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::KPLUS) {
          foundK=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
      }
      if (!foundPi && !foundK) vetoEvent;
      if (cTheta>0.6) vetoEvent;
      _cPiK->fill();
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      double sigma = _cPiK->val()*fact;
      double error = _cPiK->err()*fact;
      Estimate1DPtr cross;
      book(cross, 1, 1, 1);
      for (auto& b : cross->bins()) {
        if (inRange(sqrtS(), b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPiK;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1994_I372230);

}
