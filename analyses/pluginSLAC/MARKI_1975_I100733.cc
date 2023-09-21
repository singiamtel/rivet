// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class MARKI_1975_I100733 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(MARKI_1975_I100733);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons");
      book(_c_muons, "/TMP/sigma_muons");
      if (isCompatibleWithSqrtS(3*GeV)) {
      	book(_h_charged, 3, 1, 1);
      }
      else if (isCompatibleWithSqrtS(4.8*GeV)) {
        book(_h_charged, 3, 1, 2);
      }
      else if (isCompatibleWithSqrtS(7.4*GeV)) {
        book(_h_charged, 3, 1, 3);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
      	_c_muons->fill(); // mu+mu- + photons
      }
      // everything else
      else {
      	_c_hadrons->fill();
        if(_h_charged) {
          for (const Particle& p : fs.particles()) {
            if(PID::isCharged(p.pid())) {
              double x = 2.*p.p3().mod()/sqrtS();
              _h_charged->fill(x);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if (_h_charged) {
      	scale(_h_charged, crossSection()/ sumOfWeights() /microbarn*sqr(sqrtS()));
      }
      // R
      const double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale({_c_hadrons, _c_muons}, fact);
      BinnedEstimatePtr<string> hadrons;
      book(hadrons, 1,1,1);
      BinnedEstimatePtr<string> muons;
      book(muons, "sigma_muons", hadrons->xEdges());
      for (auto& b : hadrons->bins()) {
        const double sqs = std::stod(b.xEdge());
        if (isCompatibleWithSqrtS(sqs)) {
      	  b.set(_c_hadrons->val(), _c_hadrons->err());
          muons->bin(b.index()).set(_c_muons->val(), _c_muons->err());
      	}
      }

      Estimate0D R = *_c_hadrons/ *_c_muons;
      BinnedEstimatePtr<string> mult;
      book(mult, 2,1,1);
      for (auto& b : mult->bins()) {
        const double sqs = std::stod(b.xEdge());
        if (isCompatibleWithSqrtS(sqs)) {
      	  b.set(R.val(), R.errPos());
      	}
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c_hadrons, _c_muons;
    Histo1DPtr _h_charged;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(MARKI_1975_I100733);


}
