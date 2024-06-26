// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class CELLO_1981_I166365 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1981_I166365);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // counters for R
      book(_c_hadrons, "/TMP/sigma_hadrons",refData(1,1,1));
      book(_c_muons,   "/TMP/sigma_muons"  ,refData(1,1,1));
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22])
        _c_muons->fill(sqrtS()/GeV);
      // everything else
      else {
        _c_hadrons->fill(sqrtS()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      divide(_c_hadrons,_c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _c_hadrons, _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CELLO_1981_I166365);


}
