// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > K0K0
  class BESIII_2021_I1866051 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1866051);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_nKSKL, 1, 1, 1);
      for (const string& en : _nKSKL.binning().edges<0>()) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }

      if(ntotal==2 && nCount[130]==1 && nCount[310]==1)
	_nKSKL->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nKSKL, crossSection()/ sumOfWeights() /picobarn);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _nKSKL;
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1866051);

}
