// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  class CLEO_1991_I29927 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEO_1991_I29927);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_B    , 1, 1, 1);
      book(_c_Bstar, 2, 1, 1);
      for (const string& en : _c_B.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      unsigned int nBstar(0);
      // Get Bottom hadrons
      const Particles bhads = select(ufs.particles(), isBottomHadron);
      // find the Bstars
      for (const Particle& p : bhads) {
        if(abs(p.pid())==513 || abs(p.pid())==523) {
          if(!p.hasDescendantWith(Cuts::pid == p.pid())) ++nBstar;
        }
      }
      if(!bhads.empty())
        _c_B->fill(_ecms);
      if(nBstar!=0)
        _c_Bstar->fill(_ecms,nBstar);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      scale(_c_B,fact);
      scale(_c_Bstar,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_B, _c_Bstar;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEO_1991_I29927);


}
