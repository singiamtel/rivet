// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class HRS_1990_I280958 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(HRS_1990_I280958);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(Beam(), "Beams");
      declare(UnstableParticles(), "UFS");
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      declare(Thrust(cfs), "Thrust");

      // Book histograms
      book(_h["X"], 3, 1, 1);
      book(_h["rap_all"], 4, 1, 1);
      book(_h["rap_charm"], 5, 1, 1);
      book(_h["rap_light"], 6, 1, 1);
      book(_wLight,"TMP/wLight");
      book(_wCharm,"TMP/wCharm");

      _axes["X"] = YODA::Axis<double>({0.0, 0.1, 0.14, 0.18, 0.22, 0.26, 0.30,
                                       0.34, 0.38, 0.42, 0.46, 0.5, 0.6, 0.7});
      _axes["rap_all"] = YODA::Axis<double>({0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75,
                                             2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75});
      _axes["rap_charm"] = YODA::Axis<double>({-3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5,
                                               0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5});
      _axes["rap_light"] = _axes["rap_charm"];

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
      if (nch<5) vetoEvent;
      // Get beams and average beam momentum
      const ParticlePair& beams = apply<Beam>(event, "Beams").beams();
      const double meanBeamMom = ( beams.first.p3().mod() +
                                   beams.second.p3().mod() ) / 2.0;
      MSG_DEBUG("Avg beam momentum = " << meanBeamMom);
      // get the thrust axes
      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      const Vector3 & axis = thrust.thrustAxis();
      // unstable particles
      const UnstableParticles& ufs = apply<UnstableParticles>(event, "UFS");
      Particle pTag;
      // get the tags
      Particles Dstar = ufs.particles(Cuts::abspid==413);
      bool charmTagged = !Dstar.empty();
      if (charmTagged) {
        pTag = Dstar[0];
        for (const Particle & p : Dstar) {
          if (p.E()>pTag.E()) pTag=p;
        }
        _wCharm->fill();
      }
      bool lightTagged = false;
      if (!charmTagged) {
        for (const Particle& p : cfs.particles()) {
          if (p.p3().mod()>9.43*GeV) {
            pTag=p;
            lightTagged=true;
            _wLight->fill();
            break;
          }
        }
      }
      // sign of hemispheres if tagged
      double sign=1.;
      if(charmTagged || lightTagged) {
        if(dot(axis,pTag.p3())<0.) sign=-1.;
      }
      // now loop over the kaons
      for (const Particle & p : ufs.particles(Cuts::pid==130 || Cuts::pid==310)) {
         double xE = p.E()/meanBeamMom;
         const double energy = p.E();
         const double momT = dot(axis, p.p3());
         discfill("X", xE);
         double rap = 0.5 * std::log((energy + momT) / (energy - momT));
         discfill("rap_all", fabs(rap));
         rap *= sign;
         if (charmTagged) {
           discfill("rap_charm", rap);
         }
         else if (lightTagged) {
           discfill("rap_light", rap);
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
      scale(_h["X"], crossSection()*sqr(sqrtS())/nanobarn/sumOfWeights());
      scale(_h["rap_all"], 1./sumOfWeights());
      scale(_h["rap_light"], 1./ *_wLight);
      scale(_h["rap_charm"], 1./ *_wCharm);
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _wLight, _wCharm;
    map<string, BinnedHistoPtr<string>> _h;
    map<string, YODA::Axis<double>> _axes;
    map<string, vector<string>> _edges;
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(HRS_1990_I280958);


}
