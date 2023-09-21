// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-
  class VENUS_1995_I392360 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(VENUS_1995_I392360);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (!inRange(sqrtS(), 1*GeV, 1.5*GeV)) {
        throw Error("Invalid CMS energy for VENUS_1995_I392360");
      }
      // bin for the angle plots
      int ibin = (sqrtS()-1.0)/0.05 + 2;
      book(_h_cTheta, ibin, 1, 1);
      book(_cPi, "/TMP/nPi");
      _cmax = ibin>2 ? 0.6 : 0.4;
      if (ibin == 2)  _axis = YODA::Axis<double>(6, 0.0, 0.6);
      else            _axis = YODA::Axis<double>(4, 0.0, 0.4);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_cTheta->xEdges();
      Particles part = apply<FinalState>(event, "FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PIPLUS) {
          foundP = true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::PIMINUS)
          foundM = true;
      }
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=_cmax)    _cPi->fill();
      if (_h_cTheta )  _h_cTheta->fill(map2string(cTheta));
    }

    string map2string(const double value) const {
      const size_t idx = _axis.index(value) - 1;
      if (idx < _edges.size())  return _edges[idx];
      return "OTHER";
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) scale(_h_cTheta, fact);
      scale(_cPi, fact);
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      for (auto& b : mult->bins()) {
        const double Ecm = std::stod(b.xEdge());
        if (isCompatibleWithSqrtS(sqrtS()/GeV, Ecm)) {
          b.set(_cPi->val(), _cPi->err());
        }
      }
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_cTheta;
    CounterPtr _cPi;
    double _cmax;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(VENUS_1995_I392360);

}
