// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class CLEO_1994_I358510 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1994_I358510);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if(inRange(sqrtS()/GeV,2.,3.25)) {
        if (sqrtS()<=2.5) {
          book(_h_cTheta,2,1,1);
        }
        else if (sqrtS()<=3.) {
          book(_h_cTheta,3,1,1);
        }
       	book(_cProton, "TMP/sigma", refData(1, 1, 1));
      }
      else {
        throw Error("Invalid CMS energy for CLEO_1994_I358510");
      }
      _axis = YODA::Axis<double>({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_h_cTheta && _edges.empty())  _edges = _h_cTheta->xEdges();
      Particles part = apply<FinalState>(event,"FS").particles();
      if(part.size()!=2) vetoEvent;
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
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=0.6)    _cProton->fill(sqrtS()/GeV);
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
      Estimate1DPtr tmp;
      book(tmp,1,1,1);
      barchart(_cProton,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_cTheta;
    Histo1DPtr _cProton;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1994_I358510);

}
