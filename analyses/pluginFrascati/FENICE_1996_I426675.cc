// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> hadrons
  class FENICE_1996_I426675 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FENICE_1996_I426675);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, 1, 1, 1);
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        const double end = sqrt(std::stod(en)*GeV);
        if (isCompatibleWithSqrtS(end)) {
          _s = en;
          break;
        }
      }
      if (_s.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	vetoEvent;
      _c_hadrons->fill(_s);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c_hadrons,crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons;
    string _s;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FENICE_1996_I426675);


}
