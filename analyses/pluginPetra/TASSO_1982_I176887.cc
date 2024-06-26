// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class TASSO_1982_I176887 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(TASSO_1982_I176887);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
        book(_c_hadrons[ix], "/TMP/sigma_hadrons_"+toString(2*ix), refData<YODA::BinnedEstimate<string>>(1+2*ix,1,1));
        book(_c_muons[ix],   "/TMP/sigma_muons_"  +toString(2*ix), refData<YODA::BinnedEstimate<string>>(1+2*ix,1,1));
        for (const string& en : _c_hadrons[ix].binning().edges<0>()) {
          const size_t idx = en.find("-");
          if(idx!=std::string::npos) {
            const double emin = std::stod(en.substr(0,idx));
            const double emax = std::stod(en.substr(idx+1,string::npos));
            if(inRange(sqrtS()/GeV, emin, emax)) {
              _ecms[ix] = en;
              break;
            }
          }
          else {
            const double end = std::stod(en)*GeV;
            if (isCompatibleWithSqrtS(end)) {
              _ecms[ix] = en;
              break;
            }
          }
        }
      }
      book(_c_hadronsI, "/TMP/sigma_hadrons_1", refData<YODA::BinnedEstimate<int>>(2,1,1));
      book(_c_muonsI,   "/TMP/sigma_muons_1"  , refData<YODA::BinnedEstimate<int>>(2,1,1));
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
	for(unsigned int ix=0;ix<2;++ix) _c_muons[ix]->fill(_ecms[ix]);
        _c_muonsI->fill(round(sqrtS()/GeV));
      }
      // everything else
      else {
	for(unsigned int ix=0;ix<2;++ix) _c_hadrons[ix]->fill(_ecms[ix]);
        _c_hadronsI->fill(round(sqrtS()/GeV));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        BinnedEstimatePtr<string> mult;
        book(mult, 1+2*ix, 1, 1);
        divide(_c_hadrons[ix],_c_muons[ix],mult);
      }
      BinnedEstimatePtr<int> multI;
      book(multI, 2, 1, 1);
      divide(_c_hadronsI,_c_muonsI,multI);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons[2], _c_muons[2];
    string _ecms[2];
    BinnedHistoPtr<int> _c_hadronsI, _c_muonsI;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(TASSO_1982_I176887);


}
