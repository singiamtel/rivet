// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Hemispheres.hh"

namespace Rivet {


  /// @brief Thrust, heavy jet mass, and y3 at 58 GeV
  class TOPAZ_1993_I361661 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TOPAZ_1993_I361661);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      declare(fs, "FS");
      declare(FastJets(fs, JetAlg::DURHAM, 0.7), "DurhamJets");
      const Thrust thrust(fs);
      declare(thrust, "Thrust");
      declare(Hemispheres(thrust), "Hemispheres");

      // Book histograms
      book(_h["thrust"], 1, 1, 1);
      book(_h["rho"]   , 2, 1, 1);
      book(_h["y23"]   , 3, 1, 1);

      _axes["thrust"] = YODA::Axis<double>(22, 0.8, 5.2);
      _axes["rho"] = YODA::Axis<double>(21, 1.0, 5.2);
      _axes["y23"] = YODA::Axis<double>(18, 1.2, 8.4);

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
      // thrust
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      fillhist("thrust", -log(1.-thrust.thrust()));
      // jet mass
      const Hemispheres& hemi = apply<Hemispheres>(event, "Hemispheres");
      fillhist("rho", -log(hemi.scaledM2high()));
      // Jets
      const FastJets& durjet = apply<FastJets>(event, "DurhamJets");
      if (numParticles>=3) {
        if (durjet.clusterSeq()) fillhist("y23", -log(durjet.clusterSeq()->exclusive_ymerge_max(2)));
      }
    }

    void fillhist(const string& label, const double value) {
      string edge = "OTHER";
      const size_t idx = _axes[label].index(value);
      if (idx && idx <= _edges[label].size()) {
        edge = _edges[label][idx-1];
      }
      _h[label]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h);
      for( auto & hist : _h) {
        for(auto & b: hist.second->bins()) {
          const size_t idx = b.index();
          b.scaleW(1./_axes[hist.first].width(idx));
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    map<string, BinnedHistoPtr<string>> _h;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(TOPAZ_1993_I361661);


}
