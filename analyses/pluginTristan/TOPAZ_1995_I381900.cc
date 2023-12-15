// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TOPAZ_1995_I381900 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TOPAZ_1995_I381900);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");

      // Book histograms
      book(_h["charged"], 1, 1, 1);
      book(_h["pi"]     , 2, 1, 1);
      book(_h["Kp"]     , 2, 1, 2);
      book(_h["proton"] , 2, 1, 3);
      book(_h["K0"]     , 3, 1, 1);
      book(_wSum,"TMP/wSum");

      _axes["charged"] = YODA::Axis<double>{22, 0.6, 5.0};
      _axes["pi"] = YODA::Axis<double>{1.17, 1.47, 1.77, 2.07, 2.345, 2.545, 2.695, 2.875, 3.195,
                                       3.46, 3.56, 3.69, 3.865, 4.05, 4.28, 4.565, 4.77, 4.89};
      _axes["Kp"] = YODA::Axis<double>{1.17, 1.47, 1.77, 2.27, 2.695, 3.09, 3.46, 3.56, 3.69, 3.865, 4.05, 4.28, 4.56};
      _axes["proton"] = YODA::Axis<double>{1.17, 1.47, 1.77, 2.27, 2.695, 3.09, 3.46, 3.56, 3.69, 3.85};
      _axes["K0"] = YODA::Axis<double>{1.575, 2.025, 2.4, 2.7, 3.0, 3.3, 3.75, 4.35};
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      if (_edges.empty()) {
        for (const auto& item : _h) {
          _edges[item.first] = item.second->xEdges();
        }
      }
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      int nch = cfs.particles().size();
      if (nch<5)  vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      _wSum->fill();

      // neutral kaons
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      for (const Particle & p : ufs.particles(Cuts::pid==130 || Cuts::pid==310)) {
        const double xi = -log(p.p3().mod()/meanBeamMom);
        fillhist("K0", xi);
      }
      // charged particles
      for (const Particle& p : cfs.particles()) {
        double xi = -log(p.p3().mod()/meanBeamMom);
        fillhist("charged", xi);
        int id = abs(p.pid());
        if (id==211) {
          fillhist("pi", xi);
        }
        else if (id==321) {
          fillhist("Kp", xi);
        }
        else if (id==2212) {
          fillhist("proton", xi);
        }
      }
    }

    void fillhist(const string& label, const double value) {
      string edge = "OTHER";
      const size_t idx = _axes[label].index(value);
      if (idx && idx <= _edges[label].size())  edge = _edges[label][idx-1];
      _h[label]->fill(edge);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h, 1./ *_wSum);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _wSum;
    map<string, BinnedHistoPtr<string>> _h;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TOPAZ_1995_I381900);

}
