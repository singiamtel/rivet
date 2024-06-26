// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief gamma gamma -> pi0pi0
  class CRYSTAL_BALL_1990_I294492 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CRYSTAL_BALL_1990_I294492);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // histos
      if (inRange(sqrtS()/GeV,0.25,1.95)) {
        book(_npipi, 1, 1, 1);
        double sMin=0.27,sMax=0.5, step=0.2;
        unsigned int ihist=1;
        while (sMin<4.1) {
          if (inRange(sqrtS()/GeV, sMin, sMax)) {
            break;
          }
          sMin=sMax;
          if (fuzzyEquals(1.3, sMin)) step=0.4;
          sMax+=step;
          ihist+=1;
        }
        if (ihist<7) book(_h_cTheta[0], 2, 1, ihist);
        if (inRange(sqrtS()/GeV,1.1,1.5)) book(_h_cTheta[1], 2, 1, 7);
      }
      else {
       	throw Error("Invalid CMS energy for CRYSTAL_BALL_1990_I294492");
      }
      _comEdge = "OTHER"s;
      for (const string& edge : _npipi.binning().edges<0>()) {
        if (isCompatibleWithSqrtS(std::stod(edge))) {
          _comEdge = edge;
          break;
        }
      }
      _axis = YODA::Axis<double>({0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_h_cTheta[0] && _edges.empty())  _edges = _h_cTheta[0]->xEdges();
      if (_h_cTheta[1] && _edges.empty())  _edges = _h_cTheta[1]->xEdges();
      Particles part = apply<FinalState>(event,"FS").particles();
      if (part.size()!=2) vetoEvent;
      for(const Particle& p : part) {
        if (p.pid()!=PID::PI0) vetoEvent;
      }
      const double cTheta = abs(part[0].momentum().z()/part[0].momentum().p3().mod());
      if (cTheta<=0.8) _npipi->fill(_comEdge);
      if (_h_cTheta[0]) _h_cTheta[0]->fill(disc(cTheta));
      if (_h_cTheta[1]) _h_cTheta[1]->fill(disc(cTheta));
    }


    string disc(const double value) const {
      size_t idx = _axis.index(value);
      if (0 < idx && idx <= _edges.size())  return _edges[idx-1];
      return "OTHER"s;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      for (unsigned int ix=0; ix<2; ++ix) {
        if (_h_cTheta[ix]) {
          scale(_h_cTheta[ix], fact);
          for(auto & b : _h_cTheta[ix]->bins()) {
            const size_t idx = b.index();
            b.scaleW(1./_axis.width(idx));
          }
        }
      }
      scale(_npipi, fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npipi, _h_cTheta[2];
    YODA::Axis<double> _axis;
    vector<string> _edges;
    string _comEdge;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CRYSTAL_BALL_1990_I294492);

}
