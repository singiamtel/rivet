// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class JADE_1998_I447560 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(JADE_1998_I447560);


    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::pT > 0.1*GeV);
      declare(cfs, "CFS");
      declare(FastJets(cfs, JetAlg::DURHAM, 0.7), "DurhamJets");

      // Thrust
      const Thrust thrust(cfs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Histos
      int offset = 0;
      switch ( int(sqrtS()/GeV) ) {

        case 44:
          offset = 0;
          book(_h["thrust"], 2+offset, 1, 1);
          book(_h["MH"],     3+offset, 1, 1);
          book(_h["BT"],     4+offset, 1, 1);
          book(_h["BW"],     5+offset, 1, 1);
          book(_h["y23"], 10, 1, 1);
          break;
        case 35:
          offset = 4;
          book(_h["thrust"], 2+offset, 1, 1);
          book(_h["MH"],     3+offset, 1, 1);
          book(_h["BT"],     4+offset, 1, 1);
          book(_h["BW"],     5+offset, 1, 1);
          book(_h["y23"], 11, 1, 1);
          break;
        case 22:
          book(_h["y23"], 12, 1, 1);
          break;
      }

      _axes["thrust"] = YODA::Axis<double>{0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16,
                                           0.18, 0.2, 0.23, 0.26, 0.3, 0.35, 0.4, 0.5};
      _axes["MH"] = YODA::Axis<double>{0.0, 0.06, 0.1, 0.14, 0.18, 0.22, 0.26, 0.3,
                                       0.34, 0.38, 0.42, 0.46, 0.5, 0.55, 0.6};
      _axes["BT"] = YODA::Axis<double>{0.0, 0.03, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16,
                                       0.18, 0.2, 0.22, 0.24, 0.27, 0.3, 0.35, 0.4};
      _axes["BW"] = YODA::Axis<double>{0.0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12,
                                       0.14, 0.16, 0.18, 0.2, 0.23, 0.26, 0.3};
      _axes["y23"] = YODA::Axis<double>{0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.04,
                                        0.05, 0.06, 0.08, 0.1, 0.13, 0.16, 0.2, 0.25, 0.4};
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
        }
      }

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // JADE hadronic event selection
      if (cfs.particles().size() < 3 ) vetoEvent;

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const Vector3 & thrustAxis = thrust.thrustAxis ();
      double theta = thrustAxis.theta();
      if ( fabs(cos(theta)) >= 0.8 ) {
        MSG_DEBUG("Failed thrust angle cut: " << fabs(cos(theta)));
        vetoEvent;
      }
      /// @todo Evis, pmiss, pbal

      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");

      const double y23 = durjet.clusterSeq()->exclusive_ymerge_max(2);

      // Make sure we don't run into a segfault by trying to fill non-existing histos
      const int s = int(sqrtS()/GeV);
      if (s == 44 || s == 35) {
        discfill("thrust", 1. - thrust.thrust());
        discfill("MH", sqrt(hemi.scaledM2high()));
        discfill("BT", hemi.Bsum());
        discfill("BW", hemi.Bmax());
      }
      discfill("y23", y23);
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      for (auto& item : _h) {
        for (auto& b : item.second->bins()) {
          b.scaleW(1./_axes[item.first].width(b.index()));
        }
      }
    }

    void discfill(const string& name, const double value) {
      string edge = "OTHER";
      const size_t idx = _axes[name].index(value);
      if (idx && idx <= _axes[name].numBins())  edge = _edges[name][idx-1];
      _h[name]->fill(edge);
    }


  private:

    map<string, BinnedHistoPtr<string>> _h;
    map<string, vector<string>> _edges;
    map<string, YODA::Axis<double>> _axes;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(JADE_1998_I447560, JADE_1998_S3612880);

}
