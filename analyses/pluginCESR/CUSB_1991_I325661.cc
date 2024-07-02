// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > B* X
  class CUSB_1991_I325661 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CUSB_1991_I325661);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_Bstar, 1, 1, 1);
      for (const string& en : _c_Bstar.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms = en;
            break;
          }
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      unsigned int nBstar(0);
      // Get Bottom hadrons
      const Particles bhads = ufs.particles(Cuts::abspid==513 or Cuts::abspid==523);
      // find the Bstars
      for (const Particle& p : bhads) {
        if(abs(p.pid())==513 || abs(p.pid())==523) {
          if(!p.hasDescendantWith(Cuts::pid == p.pid())) ++nBstar;
        }
      }
      if(nBstar!=0)
        _c_Bstar->fill(_ecms,nBstar);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c_Bstar,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_Bstar;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CUSB_1991_I325661);


}
