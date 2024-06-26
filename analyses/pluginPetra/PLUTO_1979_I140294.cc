// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class PLUTO_1979_I140294 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1979_I140294);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons1, "/TMP/sigma_hadrons1", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons1,   "/TMP/sigma_muons1"  , refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_hadrons2, "/TMP/sigma_hadrons2", refData<YODA::BinnedEstimate<int>>(2,1,1));
      book(_c_muons2,   "/TMP/sigma_muons2"  , refData<YODA::BinnedEstimate<int>>(2,1,1));
      for (const string& en : _c_hadrons1.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
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
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22]) {
	_c_muons1->fill(_ecms);
	_c_muons2->fill(round(sqrtS()/GeV));
      }
      // everything else
      else {
	_c_hadrons1->fill(_ecms);
	_c_hadrons2->fill(round(sqrtS()/GeV));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> mult1;
      book(mult1, 1, 1, 1);
      divide(_c_hadrons1,_c_muons1,mult1);
      BinnedEstimatePtr<int> mult2;
      book(mult2, 2, 1, 1);
      divide(_c_hadrons2,_c_muons2,mult2);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons1, _c_muons1;
    BinnedHistoPtr<int>    _c_hadrons2, _c_muons2;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1979_I140294);


}
