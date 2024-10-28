// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  class TASSO_1990_I294755 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1990_I294755);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs(Cuts::pT >=  0.1/GeV);
      declare(cfs, "CFS");

      // Thrust and sphericity
      declare(Thrust(cfs), "Thrust");
      declare(Sphericity(cfs), "Sphericity");

      // Histos
      for (double eVal : allowedEnergies()) {

        const string en = toString(int(eVal/MeV));
        if (isCompatibleWithSqrtS(eVal, 1e-02))  _sqs = en;

        int offset = 0;
        switch (int(eVal+0.5)) {
          case 14:
            offset = 0;
            break;
          case 22:
            offset = 1;
            break;
          case 35:
            offset = 2;
            break;
          case 44:
            offset = 3;
            break;
        }

        book(_h[en+"xp0"],        2, 1, 1+offset);
        book(_h[en+"xp1"],        3, 1, 1+offset);
        book(_h[en+"xi"],         4, 1, 1+offset);
        book(_h[en+"pT"],         5, 1, 1+offset);
        book(_n[en+"sphericity"], 6, 1, 1+offset);
        book(_n[en+"aplanarity"], 7, 1, 1+offset);
        book(_n[en+"thrust"],     8, 1, 1+offset);
        book(_c[en], "/TMP/_sumWPassed_"+en);
      }
      if (_sqs == "" && !merging()) {
        throw BeamError("Invalid beam energy for " + name() + "\n");
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const ChargedFinalState& cfs = apply<ChargedFinalState>(event, "CFS");

      // TASSO hadronic event selection TODO: move this into a trigger definition
      // See page 2 in publication
      // Condition 1)  --- require at least 5 (4) 'good' tracks
      int nch = cfs.particles().size();
      if ( (int(sqrtS()/GeV) > 27 && nch < 5) || (int(sqrtS()/GeV) <= 27 && nch < 4 ) ) {
        MSG_DEBUG("Failed # good tracks cut: " << nch);
        vetoEvent;
      }
      // Condition 2) ---
      // Condition 5) --- scalar momentum (not pT!!!) sum >= 0.265*s
      double momsum = 0.0;
      for (const Particle& p : cfs.particles()) {
        const double mom = p.p3().mod();
        momsum += mom;
      }
      if (momsum <=0.265 * sqrtS()/GeV) {
        MSG_DEBUG("Failed pTsum cut: " << momsum << " < " << 0.265 * sqrtS()/GeV);
        vetoEvent;
      }

      // Raise counter for events that pass trigger conditions
      _c[_sqs]->fill();

      const Thrust& thrust = apply<Thrust>(event, "Thrust");
      //const Vector3 & thrustAxis = thrust.thrustAxis ();
      //double theta = thrustAxis.theta();
      //if ( fabs(cos(theta)) >= 0.8 ) {
        //MSG_DEBUG("Failed thrust angle cut: " << fabs(cos(theta)));
        //vetoEvent;
      //}

      const Sphericity& sphericity = apply<Sphericity>(event, "Sphericity");

      // Fill histograms in order of appearance in paper
      for (const Particle& p : cfs.particles()) {
        // Get momentum and energy of each particle.
        const Vector3 mom3 = p.p3();
        // Scaled momenta.
        const double mom = mom3.mod();
        const double scaledMom = 2.*mom/sqrtS();
        const double pTin = dot(mom3, sphericity.sphericityMajorAxis());
        const double pTout = dot(mom3, sphericity.sphericityMinorAxis());
        const double pT = sqrt(sqr(pTin)+sqr(pTout));
        _h[_sqs+"xp0"]->fill(scaledMom);
        _h[_sqs+"xp1"]->fill(scaledMom);
        _h[_sqs+"xi"]->fill(-log(scaledMom));
        _h[_sqs+"pT"]->fill(pT);
      }
      // event shapes
      _n[_sqs+"sphericity"]->fill(sphericity.sphericity());
      _n[_sqs+"aplanarity"]->fill(sphericity.aplanarity());
      _n[_sqs+"thrust"]->fill(thrust.thrust());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (auto& item : _h) {
        const string en = item.first.substr(0,5);
        if (_c[en]->sumW())  scale(item.second, 1./ *_c[en]);
      }
      normalize(_n);
    }

    /// @}


  private:

    /// @name Histograms
    /// @{
    map<string,Histo1DPtr> _h, _n;
    map<string,CounterPtr> _c;

    string _sqs = "";
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(TASSO_1990_I294755, TASSO_1990_S2148048);

}
