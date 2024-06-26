// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief R measurement
  class LENA_1982_I179431 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(LENA_1982_I179431);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for (size_t i=0; i<2; ++i) {
        const auto& ref = refData(i+1, 1, 1);
        book(_h_hadrons[i], "d0"+to_string(i+1)+"_sigma_hadrons", ref);
        book(_h_muons[i], "d0"+to_string(i+1)+"_sigma_muons", ref);
      }
      book(_c_hadrons, "d03_sigma_hadrons", refData<YODA::BinnedEstimate<string>>(3,1,1));
      book(_c_muons, "d03_sigma_muons"   , refData<YODA::BinnedEstimate<string>>(3,1,1) );
      for (const string& en : _c_hadrons.binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty() && !inRange(sqrtS()/GeV,7.4,9.43))
        MSG_ERROR("Beam energy incompatible with analysis.");
      if(_ecms.empty()) _ecms="OTHER";
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
      if (nCount[-13]==1 and nCount[13]==1 && ntotal==2+nCount[22]) {
        _h_muons[0]->fill(sqrtS()/GeV);
        _h_muons[1]->fill(sqrtS()/GeV);
        _c_muons->fill(_ecms);
      }
      // everything else
      else {
        _h_hadrons[0]->fill(sqrtS()/GeV);
        _h_hadrons[1]->fill(sqrtS()/GeV);
        _c_hadrons->fill(_ecms);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix) {
        Estimate1DPtr mult;
        book(mult, ix+1, 1, 1);
        divide(_h_hadrons[ix],_h_muons[ix],mult);
      }
      BinnedEstimatePtr<string> mult;
      book(mult, 3, 1, 1);
      divide(_c_hadrons,_c_muons,mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_hadrons, _c_muons;
    Histo1DPtr _h_hadrons[2], _h_muons[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(LENA_1982_I179431);


}
