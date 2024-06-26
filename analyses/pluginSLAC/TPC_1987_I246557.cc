// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> p pbar
  class TPC_1987_I246557 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TPC_1987_I246557);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // book histos
      if (inRange(sqrtS()/GeV,2.,2.8)) {
        if (inRange(sqrtS()/GeV,2.1,2.4)) {
          book(_h_cTheta,3,1,1);
        }
        else if(sqrtS()<=2.8) {
          book(_h_cTheta,3,1,2);
        }
        book(_cProton, 1, 1, 1);
      }
      else {
        throw Error("Invalid CMS energy for TPC_1987_I246557");
      }

      // Determine which sqrtS bin needs to be filled
      _edge = "OTHER"s;
      for (const string& edge : _cProton.binning().edges<0>()) {
        if (isCompatibleWithSqrtS(std::stod(edge)*GeV)) {
          _edge = edge;
          break;
        }
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      double cTheta(0.);
      bool foundP(false), foundM(false);
      for (const Particle& p : part) {
      if (p.pid()==PID::PROTON) {
        foundP=true;
        cTheta = abs(p.momentum().z()/p.momentum().p3().mod());
      }
      else if (p.pid()==PID::ANTIPROTON)
        foundM=true;
      }
      if (!foundP || !foundM) vetoEvent;
      if (cTheta<=0.6)    _cProton->fill(_edge);
      if (_h_cTheta) _h_cTheta->fill(cTheta);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      if (_h_cTheta) scale(_h_cTheta,fact);
      scale(_cProton, fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_cTheta;
    BinnedHistoPtr<string> _cProton;
    string _edge;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TPC_1987_I246557);

}
