// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class HRS_1992_I339573 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1992_I339573);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_h_lambda, 1 , 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty())  _edges = _h_lambda->xEdges();
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
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);

      // Final state to get particle spectra
      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==3122)) {
        double xE = p.E()/meanBeamMom;
        _h_lambda->fill(map2string(xE));
      }

    }

    string map2string(const double val) const {
      const size_t idx = _axis.index(val)-1;
      if (idx < _edges.size())  return _edges[idx];
      return "OTHER";
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double fact = sqr(sqrtS())/GeV2*crossSection()/nanobarn/sumOfWeights();
      scale(_h_lambda, fact);
      for(auto & b: _h_lambda->bins()) {
        const size_t idx = b.index();
        b.scaleW(1./_axis.width(idx));
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _h_lambda;
    YODA::Axis<double> _axis{0.0, 0.15, 0.2, 0.25, 0.3, 0.4, 0.55, 0.75};
    vector<string> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(HRS_1992_I339573);


}
