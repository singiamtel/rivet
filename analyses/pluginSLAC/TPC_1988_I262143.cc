// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TPC_1988_I262143 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TPC_1988_I262143);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_h["z_pi"],  1, 1, 1);
      book(_h["z_K"],   1, 1, 2);
      book(_h["z_p"],   1, 1, 3);
      book(_h["z_all"], 1, 1, 4);

      book(_h["z2_pi"], 5, 1, 1);
      book(_h["z2_K"],  5, 1, 2);
      book(_h["z2_p"],  5, 1, 3);

      tribook("_pi", 6, 1, 1);
      tribook("_K",  6, 1, 2);
      tribook("_p",  6, 1, 3);
      tribook("2_K", 7, 1, 1);
      tribook("2_p", 7, 1, 2);
      tribook("3_p", 7, 1, 3);
    }

    void tribook(const string& label, unsigned int d, unsigned int x, unsigned int y) {
      book(_nd["n_"+label], "TMP/n" + label, refData<YODA::BinnedEstimate<string>>(d, x, y));
      book(_nd["d_"+label], "TMP/d" + label, refData<YODA::BinnedEstimate<string>>(d, x, y));
      book(_ratio[label], d, x, y);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        for (const auto& item : _ratio) {
          _edges[item.first] = _nd["n_" + item.first]->xEdges();
        }
      }

      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = (beams.first.p3().mod() + beams.second.p3().mod()) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      for (const Particle& p : fs.particles()) {
        double xP = p.p3().mod()/meanBeamMom;
        _h["z_all"]->fill(xP);
        histfill("d_pi", xP);
        histfill("d_K",  xP);
        histfill("d_p",  xP);
        int id = abs(p.pid());
        if (id==211) {
          _h["z_pi"]->fill(xP);
          _h["z2_pi"]->fill(xP, xP);
          histfill("n_pi", xP, 100.);
          histfill("d2_K", xP);
          histfill("d2_p", xP);
          histfill("d3_p", xP);
        }
        else if (id==321) {
          _h["z_K"]->fill(xP);
          _h["z2_K"]->fill(xP, xP);
          histfill("n_K",  xP, 100.);
          histfill("n2_K", xP);
          histfill("d3_p", xP);
        }
        else if (id==2212) {
          _h["z_p"]->fill(xP);
          _h["z2_p"]->fill(xP, xP);
          histfill("n_p",  xP, 100.);
          histfill("n2_p", xP);
          histfill("n3_p", xP);
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1./sumOfWeights());
      for (auto& item : _ratio) {
        divide(_nd["n_"+item.first], _nd["d_"+item.first], item.second);
      }
    }

    /// @}

    void histfill(const string& label, const double value, const double weight = 1.0) {
      string edge = "OTHER";
      const string tag = label.substr(1);
      const size_t idx = _axis.index(value);
      if (tag == "_pi") {
        if (0.035 <= value && value <= 0.07)  edge = _edges[label][idx-3];
        if (0.09  <= value && value <= 0.11)  edge = _edges[label][idx-7];
        if (0.16  <= value && value <= 0.7 )  edge = _edges[label][idx-11];
      }
      if (tag == "_K") {
        if (0.035 <= value && value <= 0.07)  edge = _edges[label][idx-3];
        if (0.09  <= value && value <= 0.14)  edge = _edges[label][idx-7];
        if (0.25  <= value && value <= 0.7 )  edge = _edges[label][idx-12];
      }
      if (tag == "_p" || tag == "2_p") {
        if (0.035 <= value && value <= 0.11)  edge = _edges[label][idx-3];
        if (0.25  <= value && value <= 0.7 )  edge = _edges[label][idx-11];
      }
      if (tag == "2_K") {
        if (0.025 <= value && value <= 0.11)  edge = _edges[label][idx];
        if (0.25  <= value && value <= 0.7 )  edge = _edges[label][idx-8];
      }
      if (tag == "3_p") {
        if (0.035 <= value && value <= 0.07)  edge = _edges[label][idx-3];
        if (0.09  <= value && value <= 0.11)  edge = _edges[label][idx-7];
        if (0.25  <= value && value <= 0.7 )  edge = _edges[label][idx-15];
      }
      _nd[label]->fill(edge, weight);
    }


    /// @name Histograms
    /// @{
    map<string, Histo1DPtr> _h;
    map<string, BinnedHistoPtr<string>> _nd;
    map<string, BinnedEstimatePtr<string>> _ratio;
    map<string, vector<string>> _edges;
    YODA::Axis<double> _axis{0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.065,
                             0.07, 0.075, 0.08, 0.085, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14,
                             0.16, 0.18, 0.2, 0.22, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6, 0.7};


    /// }


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(TPC_1988_I262143);


}
