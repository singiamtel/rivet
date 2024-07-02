// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class ARGUS_1989_I267759 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1989_I267759);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,2.,2.9)) {
        book(_h_cTheta, 1, 1, 1);
       	book(_cProton,"TMP/nProton");
      }
      else {
        throw Error("Invalid CMS energy for ARGUS_1989_I267759");
      }
      _axis = YODA::Axis<double>({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_cTheta->xEdges();
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid() == PID::PROTON) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if(p.pid()==PID::ANTIPROTON) {
          foundM=true;
        }
      }
      if (!foundP || !foundM)  vetoEvent;
      if (cTheta<=0.6)    _cProton->fill();
      if (_h_cTheta) _h_cTheta->fill(disc(cTheta));
    }


    string disc(const double value) const {
      size_t idx = _axis.index(value);
      if (0 < idx && idx <= _axis.numBins())  return _edges[idx-1];
      return "OTHER"s;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta) {
        scale(_h_cTheta,fact);
        for(auto & b : _h_cTheta->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axis.width(idx));
        }
      }
      scale(_cProton, fact);
      Estimate1DPtr mult;
      book(mult, 2, 1, 1);
      for (auto& b : mult->bins()) {
        if (inRange(sqrtS(), b.xMin(), b.xMax())) {
          b.setVal(_cProton->val());
          b.setErr(_cProton->err());
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_cTheta;
    CounterPtr _cProton;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1989_I267759);

}
