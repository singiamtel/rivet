// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Tools/ParticleName.hh"

namespace Rivet {


  /// Measurement of charged pion, kaon, and proton production in 13 TeV pp collisions
  class CMS_2017_I1608166 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMS_2017_I1608166);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::absrap < 1.);
      declare(cfs, "CFS");
      //
      // pt spectra
      book(_h[PID::PIPLUS],  "d01-x01-y01");
      book(_h[PID::KPLUS],   "d01-x01-y02");
      book(_h[PID::PROTON],  "d01-x01-y03");
      book(_h[PID::PIMINUS], "d02-x01-y01");
      book(_h[PID::KMINUS],  "d02-x01-y02");
      book(_h[PID::PBAR],    "d02-x01-y03");

      // negative/positive ratios
      book(_s["pi-/pi+"], "d43-x01-y01");
      book(_s["k-/k+"],   "d44-x01-y01");
      book(_s["p~/p"],    "d45-x01-y01");

      // k/pi and p/pi ratios
      book(_hkpi[PID::PIPLUS], "TMP/hkpi/pi", refData<YODA::BinnedEstimate<string>>(46, 1, 1));
      book(_hkpi[PID::KPLUS],  "TMP/hkpi/k",  refData<YODA::BinnedEstimate<string>>(46, 1, 1));
      book(_hppi[PID::PIPLUS], "TMP/hppi/pi", refData<YODA::BinnedEstimate<string>>(47, 1, 1));
      book(_hppi[PID::PROTON], "TMP/hppi/p",  refData<YODA::BinnedEstimate<string>>(47, 1, 1));
      book(_s["k/pi"],    "d46-x01-y01");
      book(_s["p/pi"],    "d47-x01-y01");

      _axes[PID::PIPLUS] = YODA::Axis<double>(22, 0.1, 1.2);
      _axes[PID::KPLUS]  = YODA::Axis<double>(17, 0.2, 1.05);
      _axes[PID::PROTON] = YODA::Axis<double>(26, 0.4, 1.7);
      _axes[4] = YODA::Axis<double>(17, 0.2, 1.05);
      _axes[5] = YODA::Axis<double>(16, 0.4, 1.2);
    }


    void analyze(const Event& event) {

      if (_edges[PID::PIPLUS].empty())  _edges[PID::PIPLUS] = _h[PID::PIPLUS]->xEdges();
      if (_edges[PID::KPLUS].empty())   _edges[PID::KPLUS]  = _h[PID::KPLUS]->xEdges();
      if (_edges[PID::PROTON].empty())  _edges[PID::PROTON] = _h[PID::PROTON]->xEdges();
      if (_edges[4].empty())  _edges[4] = _hkpi[PID::PIPLUS]->xEdges();
      if (_edges[5].empty())  _edges[5] = _hppi[PID::PIPLUS]->xEdges();

      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {

        // protections against MC generators decaying long-lived particles
        if (p.hasAncestorWith(Cuts::pid == 310)  || p.hasAncestorWith(Cuts::pid == -310)  ||  // K0s
            p.hasAncestorWith(Cuts::pid == 130)  || p.hasAncestorWith(Cuts::pid == -130)  ||  // K0l
            p.hasAncestorWith(Cuts::pid == 3322) || p.hasAncestorWith(Cuts::pid == -3322) ||  // Xi0
            p.hasAncestorWith(Cuts::pid == 3122) || p.hasAncestorWith(Cuts::pid == -3122) ||  // Lambda
            p.hasAncestorWith(Cuts::pid == 3222) || p.hasAncestorWith(Cuts::pid == -3222) ||  // Sigma+/-
            p.hasAncestorWith(Cuts::pid == 3312) || p.hasAncestorWith(Cuts::pid == -3312) ||  // Xi-/+
            p.hasAncestorWith(Cuts::pid == 3334) || p.hasAncestorWith(Cuts::pid == -3334))    // Omega-/+
          continue;

        if (theParticles.find(p.pid()) != theParticles.end()) {
          // fill pt spectra
          _h[p.pid()]->fill(map2string(p.pt() / GeV, p.abspid()));
          // fill tmp histos for ratios
          if (p.abspid() != PID::PROTON) {
            _hkpi[p.abspid()]->fill(map2string(p.pt() / GeV, 4, true));
          }
          if (p.abspid() != PID::KPLUS) {
            _hppi[p.abspid()]->fill(map2string(p.pt() / GeV, 4, true));
          }
        }
      }

    }

    string map2string(const double value, const int type, const bool isRatio = false) const {
      int id = type;
      if (isRatio && id == 4)  id = PID::KPLUS;
      else if (isRatio && id == 5)  id = PID::PROTON;
      const size_t idx = _axes.at(id).index(value);
      if (idx && idx <= _edges.at(id).size()) {
        return _edges.at(id)[idx-1];
      }
      return "OTHER";
    }


    void finalize() {

      divide(_h[PID::PIMINUS], _h[PID::PIPLUS], _s["pi-/pi+"]);
      divide(_h[PID::KMINUS],  _h[PID::KPLUS],  _s["k-/k+"]);
      divide(_h[PID::PBAR],    _h[PID::PROTON], _s["p~/p"]);

      divide(_hkpi[PID::KPLUS],  _hkpi[PID::PIPLUS], _s["k/pi"]);
      divide(_hppi[PID::PROTON], _hppi[PID::PIPLUS], _s["p/pi"]);

      scale(_h, 1./2./sumOfWeights());

      for (auto& item : _h) {
        const auto& axis = _axes.at( fabs(item.first) );
        for (auto& b : item.second->bins()) {
          b.scaleW( 1.0/axis.width(b.index()) );
        }
      }

    }


    set<int> theParticles = {PID::PIPLUS, PID::KPLUS, PID::PROTON, PID::PIMINUS, PID::KMINUS, PID::PBAR};

    map<int, BinnedHistoPtr<string>> _h;
    map<int, BinnedHistoPtr<string>> _hkpi, _hppi;
    map<string, BinnedEstimatePtr<string>> _s;
    map<int, YODA::Axis<double>> _axes;
    map<int, vector<string>> _edges;

  };


  RIVET_DECLARE_PLUGIN(CMS_2017_I1608166);

}
