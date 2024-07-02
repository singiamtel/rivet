// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- > KS0 KL0 pi0
  class SND_2018_I1637194 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2018_I1637194);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKKpi, 1, 1, 1);
      for (const string& en : _nKKpi.binning().edges<0>()) {
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
      if (_ecms.empty())
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

      if(ntotal==3 && nCount[130]==1 && nCount[310]==1 && nCount[111] ==1)
        _nKKpi->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nKKpi, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nKKpi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2018_I1637194);


}
