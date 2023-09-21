// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi+pi-
  class CELLO_1992_I345437 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1992_I345437);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Final state
      declare(FinalState(),"FS");
      // check CMS energy in range
      if (sqrtS() < 0.75*GeV || sqrtS() > 2*GeV)
        throw Error("Invalid CMS energy for CELLO_1992_I345437");
      int ibin = (sqrtS()-0.70)/0.05;
      if (ibin > 0 && ibin < 19)
        book(_h_cTheta,2,1,ibin);
      if (inRange(sqrtS()/GeV,0.85,.95))
        book(_h_cTheta2,2,1,19);
      else if (inRange(sqrtS()/GeV,1.15,1.25))
        book(_h_cTheta2,2,1,20);
      else if (inRange(sqrtS()/GeV,1.25,1.35))
        book(_h_cTheta2,2,1,21);

      book(_cPi, "/TMP/nPi");
      thetaAxis = YODA::Axis<double>(16, 0.0, 0.8);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (ecmEdge == "")  ecmEdge = _cPi->bin(ecmAxis.index(sqrtS()/GeV)).xEdge();
      if (_h_cTheta  && thetaEdges[0].empty())  thetaEdges[0] = _h_cTheta->xEdges();
      if (_h_cTheta2 && thetaEdges[1].empty())  thetaEdges[1] = _h_cTheta2->xEdges();

      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size() != 2) vetoEvent;
      double cTheta(0.);
      bool foundP(false), foundM(false);
      for (const Particle& p : part) {
        if (p.pid()==PID::PIPLUS) {
          foundP=true;
          cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
        }
        else if (p.pid()==PID::PIMINUS) {
          foundM=true;
        }
      }
      if (!foundP || !foundM)  vetoEvent;
      if (cTheta <= 0.6)  _cPi->fill(ecmEdge);
      if (_h_cTheta )  _h_cTheta ->fill(map2string(cTheta, 0));
      if (_h_cTheta2)  _h_cTheta2->fill(map2string(cTheta, 1));
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta )  scale(_h_cTheta, fact);
      if (_h_cTheta2)  scale(_h_cTheta2, fact);
      scale(_cPi, fact);
    }

    /// @}

    string map2string(const double theta, const size_t alt) const {
      const size_t idx = thetaAxis.index(theta) - 1;
      if (idx < thetaEdges[alt].size())  return thetaEdges[alt][idx];
      return "OTHER";
    }

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _cPi, _h_cTheta, _h_cTheta2;
    vector<string> thetaEdges[2];
    YODA::Axis<double> thetaAxis;
    YODA::Axis<double> ecmAxis{0.75, 0.775, 0.8, 0.825, 0.85, 0.875, 0.9, 0.925, 0.95, 0.975, 1.0,
                               1.025, 1.05, 1.075, 1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3,
                               1.325, 1.35, 1.375, 1.4, 1.425, 1.45, 1.475, 1.51875, 1.6, 1.7, 1.8, 1.9, 2.0};
    string ecmEdge = "";
    /// @}

  };


  RIVET_DECLARE_PLUGIN(CELLO_1992_I345437);

}
