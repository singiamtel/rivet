// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class PLUTO_1982_I166799 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(PLUTO_1982_I166799);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      // counters for R
      book(_c_hadrons[0], "/TMP/sigma_hadrons", refData<YODA::BinnedEstimate<string>>(1,1,1));
      book(_c_hadrons[1], 2,1,1);
      book(_c_muons, "/TMP/sigma_muons",refData<YODA::BinnedEstimate<string>>(1,1,1));
      for(unsigned int ix=0;ix<2;++ix) {
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
      if(_ecms[0].empty() && _ecms[1].empty())
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
      if(nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        // mu+mu- + photons
        _c_muons->fill(_ecms[0]);
      }
      else {
        // everything else
        _c_hadrons[0]->fill(_ecms[0]);
        _c_hadrons[1]->fill(_ecms[1]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      scale(_c_hadrons[1],fact);
      BinnedEstimatePtr<string> mult;
      book(mult, 1, 1, 1);
      divide(_c_hadrons[0],_c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons[2], _c_muons;
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(PLUTO_1982_I166799);


}
