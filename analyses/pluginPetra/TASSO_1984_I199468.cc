// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TASSO_1984_I199468 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1984_I199468);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons1, "/TMP/sigma_hadrons1", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons1,   "/TMP/sigma_muons1"  , refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_hadrons2, "/TMP/sigma_hadrons2", refData(2,1,1));
      book(_c_muons2,   "/TMP/sigma_muons2"  , refData(2,1,1));
      for (const string& en : _c_hadrons1.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty() &&
         !inRange(sqrtS()/GeV,39.8, 43.1) &&
         !inRange(sqrtS()/GeV,43.2, 45.2))
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
	 ntotal==2+nCount[22]) {
	_c_muons1->fill(_ecms);
	_c_muons2->fill(sqrtS()/GeV);
      }
      // everything else
      else {
	_c_hadrons1->fill(_ecms);
	_c_hadrons2->fill(sqrtS()/GeV);
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> mult1;
      book(mult1, 1, 1, 1);
      divide(_c_hadrons1,_c_muons1,mult1);
      Estimate1DPtr mult2;
      book(mult2, 2, 1, 1);
      divide(_c_hadrons2,_c_muons2,mult2);
    }
    /// @}

    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons1, _c_muons1;
    Histo1DPtr _c_hadrons2, _c_muons2;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(TASSO_1984_I199468);


}
