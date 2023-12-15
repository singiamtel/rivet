// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Beam.hh"

namespace Rivet {


  /// @brief K+/- and protno spectra
  class DELPHI_1995_I394052 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_1995_I394052);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "FS");

      // Book histograms
      book(_h_kaon_p   , 3, 1, 1);
      book(_h_kaon_x   , 5, 1, 1);
      book(_h_proton_p , 4, 1, 1);
      book(_h_proton_x , 6, 1, 1);

      _axis[0] = YODA::Axis<double>({0.017543859649122806, 0.024122807017543862, 0.03070175438596491,
                                     0.03728070175438596,  0.043859649122807015, 0.05482456140350877,
                                     0.06578947368421052,  0.07675438596491228,  0.08771929824561403,
                                     0.21929824561403508,  0.2631578947368421,   0.30701754385964913,
                                     0.3508771929824561,   0.39473684210526316,  0.43859649122807015,
                                     0.5043859649122807});
      _axis[1] = YODA::Axis<double>({0.03070175438596491, 0.03728070175438596, 0.043859649122807015,
                                     0.05482456140350877, 0.06578947368421052, 0.07675438596491228,
                                     0.08771929824561403, 0.09868421052631579, 0.10964912280701754});
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      if (_edges[0].empty())  _edges[0] = _h_kaon_x->xEdges();
      if (_edges[1].empty())  _edges[1] = _h_proton_x->xEdges();

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
      const double meanBeamMom = ( beams.first.p3().mod() + beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      for (const Particle& p : fs.particles(Cuts::abspid==321 or Cuts::abspid==2212)) {
        double modp = p.p3().mod();
        double xp = modp/meanBeamMom;
        if(abs(p.pid())==321) {
          _h_kaon_p->fill(modp);
          _h_kaon_x->fill(map2string(xp, 0));
        }
        else {
          _h_proton_p->fill(modp);
          _h_proton_x->fill(map2string(xp, 1));
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_kaon_p  ,1./sumOfWeights());
      scale(_h_kaon_x  ,1./sumOfWeights());
      scale(_h_proton_p,1./sumOfWeights());
      scale(_h_proton_x,1./sumOfWeights());
    }

    /// @}

    string map2string(const double value, const size_t type) const {
      const size_t idx = _axis[type].index(value);
      if (idx && idx < 9)  return _edges[type][idx-1];
      if (idx > 9 && idx <= _edges[type].size()+1)  return _edges[type][idx-2];
      return "OTHER";
    }


    /// @name Histograms
    /// @{
    Histo1DPtr _h_kaon_p, _h_proton_p;
    BinnedHistoPtr<string> _h_kaon_x, _h_proton_x;
    vector<string> _edges[2];
    YODA::Axis<double> _axis[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(DELPHI_1995_I394052);


}
