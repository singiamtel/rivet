// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief Add a short analysis description here
  class CMD3_2016_I1444990 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2016_I1444990);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      declare(FinalState(), "FS");
      book(_nK0K0, "TMP/K0K0", refData(1,1,6));

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
	 nCount[130]==1 && nCount[310]==1)
	_nK0K0->fill(sqrtS()/MeV);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nK0K0,crossSection()/ sumOfWeights() /nanobarn);
      Estimate1DPtr tmp;
      book(tmp,1,1,6);
      barchart(_nK0K0,tmp);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _nK0K0;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2016_I1444990);


}
