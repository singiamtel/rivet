// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief pi0 production cross section measured at STAR
  class STAR_2014_I1253360 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(STAR_2014_I1253360);

    const double ETA_MIN = 0.8;
    const double ETA_MAX = 2.0;
    const double DELTA_ETA = ETA_MAX - ETA_MIN;

    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      const UnstableParticles fs(Cuts::eta > ETA_MIN && Cuts::eta < ETA_MAX && Cuts::pid == 111);
      declare(fs, "fs");

      // Book histograms
      book(_h_pi0_pt, 8, 1, 1);
      _axis = YODA::Axis<double>({5., 6., 7., 8., 9., 10., 12., 16.});

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty())  _edges = _h_pi0_pt->xEdges();

      const Particles& particles = apply<UnstableParticles>(event, "fs").particles();

      for (const Particle& p : particles) {
        discfill(p.pT() / GeV, 1 / (p.pT() / GeV));
      }

    }

    void discfill(const double value, double weight) {
      const size_t idx = _axis.index(value);
      string edge("OTHER");
      if (idx && idx <= _edges.size()) {
        edge = _edges[idx-1];
        weight /= _axis.width(idx);
      }
      _h_pi0_pt->fill(edge, weight);
    }

    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_pi0_pt, crossSection()/millibarn/sumW() / (2 * pi) / DELTA_ETA);

    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _h_pi0_pt;
    YODA::Axis<double> _axis;
    vector<string> _edges;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(STAR_2014_I1253360);

}
