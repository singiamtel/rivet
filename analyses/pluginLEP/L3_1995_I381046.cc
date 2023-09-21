// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief B* production
  class L3_1995_I381046 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(L3_1995_I381046);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // // Initialise and register projections
      declare(ChargedFinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_bStar, "/TMP/cbStar", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_B    , "/TMP/cB", refData<YODA::BinnedEstimate<string>>(1,1,1));

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // First, veto on leptonic events by requiring at least 4 charged FS particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      const size_t numParticles = fs.particles().size();

      // Even if we only generate hadronic events, we still need a cut on numCharged >= 2.
      if (numParticles < 2) {
        MSG_DEBUG("Failed leptonic event cut");
        vetoEvent;
      }
      MSG_DEBUG("Passed leptonic event cut");

      for (const Particle& p : apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==513 or Cuts::abspid==523 or
									       Cuts::abspid==511 or Cuts::abspid==521)) {
        // count number of Bs not from mixing or B*
        if (p.abspid()==511 || p.abspid()==521) {
          if (p.parents()[0].abspid()==p.abspid()) continue;
          if (p.parents()[0].abspid()==513 || p.parents()[0].abspid()==523) continue;
          _c_B->fill(Ecm);
        }
        else {
          _c_bStar->fill(Ecm);	// B*
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      // no of B*/B+B*
      BinnedEstimatePtr<string> h1;
      book(h1,1,1,1);
      *h1 = *_c_bStar / (*_c_bStar + *_c_B);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_bStar, _c_B;
    const string Ecm = "91.2";
    /// @}


  };


  // The hook for the plugin system
  RIVET_DECLARE_PLUGIN(L3_1995_I381046);


}
