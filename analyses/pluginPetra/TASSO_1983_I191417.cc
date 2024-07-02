// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class TASSO_1983_I191417 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1983_I191417);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,2.,3.1)) {
        if (sqrtS()<=2.4) {
          book(_h_cTheta,2,1,1);
        }
        else if (sqrtS()<=2.8) {
          book(_h_cTheta,2,1,2);
        }
       	book(_cProton, 1, 1, 1);
      }
      else {
        throw Error("Invalid CMS energy for TASSO_1983_I191417");
      }

      _ecms = "OTHER"s;
      for (const string& edge : _cProton.binning().edges<0>()) {
        if (isCompatibleWithSqrtS(std::stod(edge)*GeV)) {
          _ecms = edge;
          break;
        }
      }
      _axis = YODA::Axis<double>({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if(_h_cTheta && _edges.empty()) _edges=_h_cTheta->xEdges();
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PROTON) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::ANTIPROTON) {
          foundM=true;
        }
      }
      if (!foundP || !foundM)  vetoEvent;
      if (cTheta<=0.6) _cProton->fill(_ecms);
      if (_h_cTheta) {
        const size_t idx = _axis.index(cTheta);
        if (0 < idx && idx < _axis.numBins())
          _h_cTheta->fill(_edges[idx-1]);
        else 
          _h_cTheta->fill("OTHER"s);
      }
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
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_cTheta, _cProton;
    string _ecms;
    vector<string> _edges;
    YODA::Axis<double> _axis;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1983_I191417);

}
