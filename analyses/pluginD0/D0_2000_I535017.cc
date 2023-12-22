// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/LeptonFinder.hh"

namespace Rivet {


  /// D0 transverse momentum of the W boson
  class D0_2000_I535017 : public Analysis {
  public:

    RIVET_DEFAULT_ANALYSIS_CTOR(D0_2000_I535017);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare("MET", MissingMomentum());
      LeptonFinder ef(0.2, Cuts::abseta < 5 && Cuts::abspid == PID::ELECTRON);
      declare(ef, "Elecs");

      book(_h_W_pT ,1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const P4& pmiss = apply<MissingMom>(event, "MET").missingMom();
      const Particles& es = apply<LeptonFinder>(event, "Elecs").particles();
      const int ifound = closestMatchIndex(es, pmiss, Kin::mass, 80.4*GeV, 0*GeV, 200*GeV);

      if (ifound < 0) vetoEvent;
      _h_W_pT->fill((pmiss+es[ifound].mom()).pT()/GeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_W_pT, crossSection()/sumOfWeights());
    }

    /// @}


  private:

    /// Histogram
    Histo1DPtr _h_W_pT;

  };


  RIVET_DECLARE_ALIASED_PLUGIN(D0_2000_I535017, D0_2000_S4480767);

}
