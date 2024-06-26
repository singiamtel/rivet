// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class PLUTO_1980_I152291 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1980_I152291);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
        book(_c_hadrons[ix], "/TMP/sigma_hadrons_"+toString(ix+1), refData(1+ix,1,1));
        book(_c_muons  [ix], "/TMP/sigma_muons_"  +toString(ix+1), refData(1+ix,1,1));
      }
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
	 ntotal==2+nCount[22]) {
	for(unsigned int ix=0;ix<2;++ix) _c_muons[ix]->fill(sqrtS());
      }
      // everything else
      else {
	for(unsigned int ix=0;ix<2;++ix) _c_hadrons[ix]->fill(sqrtS());
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        Estimate1DPtr mult;
        book(mult, ix+1, 1, 1);
        divide(_c_hadrons[ix], _c_muons[ix],mult);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _c_hadrons[2], _c_muons[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1980_I152291);


}
