// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Cross-section for KSO K+-pi-+ between 3.8 and 4.6 GeV
  class BESIII_2018_I1691798 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1691798);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      book(_nKKpi, 1, 1, 1);

      for (const string& en : _nKKpi.binning().edges<0>()) {
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

      if(ntotal==3 && nCount[310]==1 &&
         ((nCount[ 321]=1 &&  nCount[-211] ==1) ||
          (nCount[-321]=1 &&  nCount[ 211] ==1)))
        _nKKpi->fill(_ecms);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nKKpi, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nKKpi;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1691798);


}
