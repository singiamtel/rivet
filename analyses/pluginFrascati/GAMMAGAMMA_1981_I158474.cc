// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > hadrons
  class GAMMAGAMMA_1981_I158474 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(GAMMAGAMMA_1981_I158474);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<6;++ix)
        book(_h[ix],"TMP/h_"+toString(ix+1),refData(1,1,1+ix));
      for(unsigned int ix=6;ix<9;++ix)
        book(_h[ix],"TMP/h_"+toString(ix+1),refData(2,1,1));
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
	_h[8]->fill(sqrtS());
      else {
	if(ntotal==3 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==1 ) {
	  _h[0]->fill(sqrtS());
	}
	if(ntotal==4 && nCount[211] == 1 && nCount[-211]==1 && nCount[111]==2 ) {
	  _h[1]->fill(sqrtS());
	}
	if(ntotal==5 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==1 ) {
	  _h[2]->fill(sqrtS());
	}
	if(ntotal==6 && nCount[211] == 2 && nCount[-211]==2 && nCount[111]==2 ) {
	  _h[3]->fill(sqrtS());
	}
	if(nCount[211] == 1 && nCount[-211]==1 && ntotal == 2+nCount[111]) {
	  _h[6]->fill(sqrtS());
	}
	if(nCount[211] == 2 && nCount[-211]==2 && ntotal == 4+nCount[111]) {
	  _h[7]->fill(sqrtS());
	}
	if((nCount[211]+nCount[-211]+nCount[111])==ntotal ) {
	  if(ntotal==3 || ntotal ==5)
	    _h[4]->fill(sqrtS());
	  else if(ntotal==4 || ntotal==6)
	    _h[5] ->fill(sqrtS());
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<9;++ix) {
        scale(_h[ix],fact);
        if(ix>=6) continue;
        Estimate1DPtr tmp;
        book(tmp,1,1,1+ix);
        barchart(_h[ix],tmp);
      }
      Estimate1DPtr tmp;
      book(tmp,2,1,1);
      divide(_h[6],_h[8],tmp);
      book(tmp,2,1,2);
      divide(_h[7],_h[8],tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[9];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(GAMMAGAMMA_1981_I158474);


}
