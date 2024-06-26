// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Measurement of R
  class PDG_R : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PDG_R);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons",refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons  , "/TMP/sigma_muons"  ,refData<YODA::BinnedEstimate<string>>(1,1,1));
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        const size_t idx1 = en.find("_");
        const string en2 = idx1!=std::string::npos ? en.substr(0,idx1) : en;
        const size_t idx2 = en.find("-");
        if(idx2!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx2));
          const double emax = std::stod(en.substr(idx2+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms.push_back(en);
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms.push_back(en);
            break;
          }
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
        for(const string & en : _ecms)
          _c_muons->fill(en);
      }
      // everything else
      else {
        for(const string & en : _ecms)
          _c_hadrons->fill(en);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      divide(_c_hadrons, _c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    vector<string> _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PDG_R);

}
