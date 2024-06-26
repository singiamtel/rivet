// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> p pbar
  class FENICE_1994_I377833 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(FENICE_1994_I377833);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nProton, 1, 1, 1);
      for (const string& en : _nProton.binning().edges<0>()) {
        const double end = sqrt(std::stod(en)*GeV);
        if (isCompatibleWithSqrtS(end)) {
          _s = en;
          break;
        }
      }
      if (_s.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
	nCount[p.pid()] += 1;
	++ntotal;
      }
      if(ntotal==2 && nCount[2212]==1 && nCount[-2212]==1)
	_nProton->fill(_s);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nProton, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nProton;
    string _s;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(FENICE_1994_I377833);


}
