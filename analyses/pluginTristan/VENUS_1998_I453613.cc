// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief charged multiplicity at 58 GeV
  class VENUS_1998_I453613 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(VENUS_1998_I453613);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");
      book(_cLight , 1, 1, 2);
      book(_cBottom, 1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      if (cfs.size() < 2) vetoEvent;


      int flavour = 0;
      const InitialQuarks& iqf = apply<InitialQuarks>(event, "IQF");

      // If we only have two quarks (qqbar), just take the flavour.
      // If we have more than two quarks, look for the highest energetic q-qbar pair.
      if (iqf.particles().size() == 2) {
        flavour = iqf.particles().front().abspid();
      }
      else {
        map<int, double> quarkmap;
        for ( const Particle& p : iqf.particles()) {
          if (quarkmap[p.pid()] < p.E()) {
            quarkmap[p.pid()] = p.E();
          }
        }
        double maxenergy = 0.;
        for (int i = 1; i <= 5; ++i) {
          if (quarkmap[i]+quarkmap[-i] > maxenergy) {
            flavour = i;
          }
        }
      }
      const size_t numParticles = cfs.particles().size();
      switch (flavour) {
      case 1: case 2: case 3:
        _cLight->fill(round(sqrtS()),numParticles);
        break;
      case 5:
        _cBottom->fill(round(sqrtS()),numParticles);
        break;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<int> hDiff;
      book(hDiff,1, 1, 3);
      if(_cBottom->bin(1).numEntries()>0 &&
         _cLight ->bin(1).numEntries()>0) {
        double val = _cBottom->bin(1).mean(2) - _cLight->bin(1).mean(2);
        double err = sqrt(sqr(_cBottom->bin(1).stdErr(2)) +
                          sqr(_cLight ->bin(1).stdErr(2)));
        hDiff->bin(1).setVal(val);
        hDiff->bin(1).setErr(err);
      }
    }

    /// @}

    /// @name Multiplicities
    /// @{
    BinnedProfilePtr<int> _cLight;
    BinnedProfilePtr<int> _cCharm;
    BinnedProfilePtr<int> _cBottom;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(VENUS_1998_I453613);


}
