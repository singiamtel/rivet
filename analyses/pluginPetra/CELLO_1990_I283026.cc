// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Strange hadrons at 35GeV
  class CELLO_1990_I283026 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1990_I283026);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h["K0"],     1, 1, 1);
      book(_h["Kstar"],  2, 1, 1);
      book(_h["Lambda"], 3, 1, 1);

      _axes["K0"] = YODA::Axis<double>({0.0315, 0.0365, 0.0535, 0.0745, 0.0955, 0.1165,
                                        0.1415, 0.175, 0.2375, 0.282, 0.3775, 0.652});
      _axes["Kstar"] = YODA::Axis<double>({0.055, 0.195, 0.335, 0.537, 0.883});
      _axes["Lambda"] = YODA::Axis<double>({0.067, 0.085, 0.103, 0.131, 0.175, 0.269, 0.447, 0.723});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
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
      const double meanBeamMom = 0.5*(beams.first.p3().mod() + beams.second.p3().mod());
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state of unstable particles to get particle spectra
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");

      for (const Particle& p : ufs.particles()) {
        const int id = p.abspid();
        if (id == PID::K0S || id == PID::K0L) {
          discfill("K0", p.E()/meanBeamMom);
        }
        else if(abs(id)==323) {
          discfill("Kstar", p.E()/meanBeamMom);
        }
        else if(abs(id)==3122) {
          discfill("Lambda", p.E()/meanBeamMom);
        }
      }
    }

    void discfill(const string& name, const double value) {
      string edge = "OTHER";
      const size_t idx = _axes[name].index(value);
      if (idx && idx <= _edges[name].size())  edge = _edges[name][idx-1];
      _h[name]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1./sumOfWeights());
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, BinnedHistoPtr<string>> _h;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CELLO_1990_I283026);


}
