// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief pi, K and proton spectra at 12 and 30 GeV
  class TASSO_1980_I153656 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1980_I153656);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      _iHist=-1;
      sqs = 1.;
      if (isCompatibleWithSqrtS(12*GeV)) {
        _iHist = 0;
        sqs = 12.;
      }
      else if (isCompatibleWithSqrtS(30*GeV)) {
        _iHist = 1;
        sqs = 30.;
      }
      else
        MSG_ERROR("Beam energy " << sqrtS() << " GeV not supported!");

      book(_h["p_pi"],3*_iHist+2,1,1);
      book(_h["x_pi"],3*_iHist+2,1,2);
      book(_h["p_K"] ,3*_iHist+3,1,1);
      book(_h["x_K"] ,3*_iHist+3,1,2);
      book(_h["p_p"] ,3*_iHist+4,1,1);
      book(_h["x_p"] ,3*_iHist+4,1,2);

      tribook("pi", 3*_iHist+ 8, 1, 1);
      tribook("K",  3*_iHist+ 9, 1, 1);
      tribook("p",  3*_iHist+10, 1, 1);

      if (_iHist) {
        _axes["pi"] = YODA::Axis<double>({0.325, 0.375, 0.425, 0.5, 0.575, 0.7, 0.9, 1.1, 1.3, 1.5});
        _axes["K"]  = YODA::Axis<double>({0.4, 0.5, 0.575, 0.7, 0.9, 1.1});
        _axes["p"]  = YODA::Axis<double>({0.5, 0.7, 1.2325, 2.0975});
        _axes["r"]  = YODA::Axis<double>({0.4, 0.5, 0.675, 0.925});
      }
      else {
        _axes["pi"] = YODA::Axis<double>({0.3, 0.4, 0.5, 0.675, 1.05, 1.55});
        _axes["K"]  = YODA::Axis<double>({0.4, 0.5, 0.675, 0.925});
        _axes["p"]  = YODA::Axis<double>({0.475, 0.725, 1.2125, 1.9375});
        _axes["r"]  = YODA::Axis<double>({0.4, 0.5, 0.575, 0.7, 0.9, 1.1});
      }
      _axes["rp"]  = YODA::Axis<double>({0.475, 0.725, 1.0, 2.2});
    }

    void tribook(const string& label, unsigned int d, unsigned int x, unsigned int y) {
      book(_h["n_"+label], "TMP/n_"+label, refData<YODA::BinnedEstimate<string>>(d, x, y));
      book(_h["d_"+label], "TMP/d_"+label, refData<YODA::BinnedEstimate<string>>(d, x, y));
      book(_r[label], d, x, y);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
        }
        for (const auto& item : _r) {
          _edges["r"+item.first] = item.second->xEdges();
        }
      }
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const ChargedFinalState& fs = apply<ChargedFinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      for (const Particle& p : fs.particles()) {
        double modp = p.p3().mod();
        fillND("d_pi", modp);
        fillND("d_K", modp);
        fillND("d_p", modp);
        double beta = modp/p.E();
        double xE = p.E()/meanBeamMom;
        if (abs(p.pid())==211) {
          fillhist("p_pi", modp);
          fillhist("x_pi", xE, 1./beta);
          fillND("n_pi", modp);
        }
        else if (abs(p.pid())==321) {
          fillhist("p_K", modp);
          fillhist("x_K", xE, 1./beta);
          fillND("n_K", modp);
        }
        else if (abs(p.pid())==2212) {
          fillhist("p_p", modp);
          fillhist("x_p", xE, 1./beta);
          fillND("n_p", modp);
        }
      }
    }

    void fillhist(const string& label, const double value, const double weight = 1.0) {
      string edge = "OTHER";
      const string tag = label.substr(2);
      const size_t idx = _axes[tag].index(value);
      if (idx && idx <= _edges[label].size())  edge = _edges[label][idx];
      _h[label]->fill(edge, weight);
    }

    void fillND(const string& label, const double value) {
      string edge = "OTHER";
      const string tag = label.substr(2);
      const size_t idx = _axes[(tag == "p")? "rp" : "r"].index(value);
      if (idx && idx <= _edges["r"+tag].size())  edge = _edges[label][idx];
      _h[label]->fill(edge);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h["p_pi"], crossSection()/nanobarn/sumOfWeights());
      scale(_h["x_pi"], sqr(sqs)*crossSection()/microbarn/sumOfWeights());
      scale(_h["p_K"],  crossSection()/nanobarn/sumOfWeights());
      scale(_h["x_K"],  sqr(sqs)*crossSection()/microbarn/sumOfWeights());
      scale(_h["p_p"],  crossSection()/nanobarn/sumOfWeights());
      scale(_h["x_p"],  sqr(sqs)*crossSection()/microbarn/sumOfWeights());

      for (auto& item : _r) {
        divide(_h["n_"+item.first], _h["d_"+item.first], item.second);
      }

    }

    /// @}


    /// @name Histograms
    /// @{
    map<string,BinnedHistoPtr<string>> _h;
    map<string,BinnedEstimatePtr<string>> _r;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    int _iHist;
    double sqs;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1980_I153656);


}
