// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R and hadron mult in e+e-
  class GAMMAGAMMA_1979_I141722 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(GAMMAGAMMA_1979_I141722);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      book(_c_hadrons, "/TMP/sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_muons  , "/TMP/sigma_muons"  , refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_charged, 2, 1, 1);
      book(_c_neutral, 2, 1, 2);
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        const size_t idx = en.find("-");
        if(idx!=std::string::npos) {
          const double emin = std::stod(en.substr(0,idx));
          const double emax = std::stod(en.substr(idx+1,string::npos));
          if(inRange(sqrtS()/GeV, emin, emax)) {
            _ecms = en;
            break;
          }
        }
        else {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms = en;
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
      int ntotal(0),ncharged(0),nneutral(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
	if(PID::isCharged(p.pid()))
	  ncharged += 1;
	else
	  nneutral += 1;
      }
      // mu+mu- + photons
      if(nCount[-13]==1 and nCount[13]==1 &&
	 ntotal==2+nCount[22])
	_c_muons->fill(_ecms);
      // everything else
      else {
	if(ntotal==2) vetoEvent;
	_c_hadrons->fill(_ecms);
	_c_charged->fill(sqrtS()/GeV,ncharged);
	_c_neutral->fill(sqrtS()/GeV,nneutral);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      BinnedEstimatePtr<string> ratio;
      book(ratio,1,1,1);
      divide(_c_hadrons,_c_muons,ratio);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    Profile1DPtr _c_neutral,_c_charged;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(GAMMAGAMMA_1979_I141722);


}
