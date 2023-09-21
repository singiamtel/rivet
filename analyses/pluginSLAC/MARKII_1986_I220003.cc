// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-/K+ K-
  class MARKII_1986_I220003 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKII_1986_I220003);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (sqrtS()<1.5*GeV || sqrtS()>5*GeV)
        throw Error("Invalid CMS energy for MARKII_1986_I220003");
      book(_cPiK[0], "/TMP/nPiK_0");
      book(_cPiK[1], "/TMP/nPiK_1");
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      if (part[0].pid()!=-part[1].pid()) vetoEvent;
      double cTheta(0.);
      bool foundPi(false),foundK(false);
      for (const Particle & p : part) {
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
      if (cTheta>0.5) vetoEvent;
      if (cTheta>0.3)
        _cPiK[1]->fill();
      else
        _cPiK[0]->fill();
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ix=0;ix<2;++ix) {
        double sigma = _cPiK[ix]->val()*fact;
        double error = _cPiK[ix]->err()*fact;
        Estimate1DPtr cross;
        book(cross, 1, 1, 1+ix);
        for (auto& b : cross->bins()) {
          if (inRange(sqrtS(), b.xMin(), b.xMax())) {
            b.set(sigma, error);
          }
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    CounterPtr _cPiK[2];
    ///@}


  };


  RIVET_DECLARE_PLUGIN(MARKII_1986_I220003);

}
