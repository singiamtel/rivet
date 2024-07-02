// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-
  class BELLE_2007_I749358 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2007_I749358);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      _energyAxis = YODA::Axis<double>(140, 0.8, 1.5);
      _thetaAxis = YODA::Axis<double>(12, 0.0, 0.6);
      if (sqrtS()<0.8*GeV || sqrtS()>1.5*GeV)
        throw Error("Invalid CMS energy for BELLE_2007_I749358");
      eIndex = _energyAxis.index(sqrtS()/GeV) - 1;

      // bin for the angle plots
      int ibin = (sqrtS()-0.8)/0.005 + 2;
      book(_h_cTheta,ibin,1,1);
      book(_cPi, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_EcmEdges.empty())  _EcmEdges = _cPi->xEdges();
      if (_thetaEdges.empty())  _thetaEdges = _h_cTheta->xEdges();

      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false),foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PIPLUS) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::PIMINUS)
          foundM=true;
      }
      if (!foundP || !foundM)  vetoEvent;
      if (cTheta<=0.6)  _cPi->fill(_EcmEdges[eIndex]);
      if (_h_cTheta ) {
        size_t idx = _thetaAxis.index(cTheta) - 1;
        _h_cTheta->fill(idx < _thetaEdges.size() ? _thetaEdges[idx] : "OTHER"s);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta ) {
        scale(_h_cTheta ,fact);
        for(auto & b : _h_cTheta->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_thetaAxis.width(idx));
        }
      }
      scale(_cPi, fact);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_cTheta, _cPi;
    YODA::Axis<double> _energyAxis;
    YODA::Axis<double> _thetaAxis;
    vector<string> _EcmEdges, _thetaEdges;
    size_t eIndex;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2007_I749358);

}
