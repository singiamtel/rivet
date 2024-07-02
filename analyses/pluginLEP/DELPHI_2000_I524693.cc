// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/Beam.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/Thrust.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/ParisiTensor.hh"
#include "Rivet/Projections/Hemispheres.hh"
#include <cmath>

#define I_KNOW_THE_INITIAL_QUARKS_PROJECTION_IS_DODGY_BUT_NEED_TO_USE_IT
#include "Rivet/Projections/InitialQuarks.hh"

namespace Rivet {


  /// @brief DELPHI multiplicities at various energies
  ///
  /// @author Peter Richardson
  class DELPHI_2000_I524693 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(DELPHI_2000_I524693);


    /// @name Analysis methods
    /// @{

    void init() {
      // Projections
      declare(Beam(), "Beams");
      declare(ChargedFinalState(), "CFS");
      declare(InitialQuarks(), "IQF");

      // Histograms
      book(_hLight,  1,1,3);
      book(_hCharm,  1,1,2);
      book(_hBottom, 1,1,1);

    }


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
        for (const Particle& p : iqf.particles()) {
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
        _hLight ->fill(int(sqrtS()),numParticles);
        break;
      case 4:
        _hCharm ->fill(int(sqrtS()),numParticles);
        break;
      case 5:
        _hBottom->fill(int(sqrtS()),numParticles);
        break;
      }
    }


    void finalize() {
      BinnedEstimatePtr<int> hDiff;
      book(hDiff,1,1,4);
      for(unsigned int ix=0;ix<hDiff->numBins();++ix) {
        if(_hBottom->bin(ix+1).numEntries()>0 &&
           _hLight ->bin(ix+1).numEntries()>0) {
          double val = _hBottom->bin(ix+1).mean(2) - _hLight->bin(ix+1).mean(2);
          double err = sqrt(sqr(_hBottom->bin(ix+1).stdErr(2)) +
                            sqr(_hLight ->bin(ix+1).stdErr(2)));
          hDiff->bin(ix+1).set(val,err);
        }
      }
    }

    /// @}


  private:

    vector<Estimate1DPtr> _mult;

    /// @name Multiplicities
    /// @{
    BinnedProfilePtr<int>  _hLight,_hCharm,_hBottom;
    /// @}

  };



  RIVET_DECLARE_ALIASED_PLUGIN(DELPHI_2000_I524693, DELPHI_2000_S4328825);

}
