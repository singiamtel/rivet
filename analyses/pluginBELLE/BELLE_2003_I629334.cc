// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> K+K-
  class BELLE_2003_I629334 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2003_I629334);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (sqrtS()<1.4*GeV || sqrtS()>2.4*GeV)
        throw Error("Invalid CMS energy for BELLE_2003_I629334");
      // bin for the angle plots
      int ibin = (sqrtS()-1.40)/0.04;
      book(_h_cTheta,2+ibin/4,1,ibin%4+1);
      book(_cK, "/TMP/nK");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::KPLUS) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::KMINUS)
          foundM=true;
      }
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=0.6)    _cK->fill();
      if (_h_cTheta ) _h_cTheta ->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta ,fact);
      double sigma = _cK->val()*fact;
      double error = _cK->err()*fact;
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS(), b.xMin(), b.xMax())) {
          b.set(sigma, error);
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    Histo1DPtr _h_cTheta;
    CounterPtr _cK;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2003_I629334);

}
