// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+ e- > KS0 KL0
  class CMD2_1999_I502164 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD2_1999_I502164);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      for(unsigned int ix=0;ix<4;++ix)
        book(_sigma[ix], "TMP/sigma_"+toString(ix), refData(1+ix, 1, 1));

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
      if(ntotal==2 &&
	 nCount[130]==1 && nCount[310]==1) {
	for(unsigned int ix=0;ix<4;++ix) _sigma[ix]->fill(sqrtS()/MeV);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<4;++ix) {
        scale(_sigma[ix],fact);
        Estimate1DPtr tmp;
        book(tmp, 1+ix, 1, 1);
        barchart(_sigma[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD2_1999_I502164);


}
