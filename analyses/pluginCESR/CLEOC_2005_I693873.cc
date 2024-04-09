// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+ e- > pi+pi-, K+K- and p pbar at 3.671 GeV
  class CLEOC_2005_I693873 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOC_2005_I693873);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_npipi , 1,1,1);
      book(_nKK   , 1,1,2);
      book(_nppbar, 1,1,3);
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
      if(ntotal!=2) vetoEvent;

      if(nCount[211]==1 && nCount[-211]==1)
	_npipi->fill(ecms);
      else if(nCount[321]==1 && nCount[-321]==1)
	_nKK->fill(ecms);
      else if(nCount[2212]==1 && nCount[-2212]==1)
	_nppbar->fill(ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_npipi ,crossSection()/ sumOfWeights() /picobarn);
      scale(_nKK   ,crossSection()/ sumOfWeights() /picobarn);
      scale(_nppbar,crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npipi,_nKK,_nppbar;
    string ecms = "3.671";
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOC_2005_I693873);


}
