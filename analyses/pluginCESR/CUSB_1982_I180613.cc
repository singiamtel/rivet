// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class CUSB_1982_I180613 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CUSB_1982_I180613);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons",refData(3,1,1));
      book(_c_muons, "/TMP/sigma_muons"    ,refData(3,1,1));
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
      if(nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _c_muons->fill(sqrtS()/GeV);
      }
      else { // everything else
        _c_hadrons->fill(sqrtS()/GeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate1DPtr mult;
      book(mult, 3, 1, 1);
      divide(_c_hadrons,_c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _c_hadrons, _c_muons;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CUSB_1982_I180613);


}
