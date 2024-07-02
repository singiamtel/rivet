// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class SND_2001_I579319 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2001_I579319);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<2;++ix) {
        book(_npion[ix], 1+ix, 1, 1);
        for (const string& en : _npion[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty())
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
      if(ntotal!=4) vetoEvent;
      if(nCount[-211]==2&&nCount[211]==2)
	_npion[0]->fill(_ecms[0]);
      else if(nCount[-211]==1&&nCount[211]==1&&nCount[111]==2)
	_npion[1]->fill(_ecms[1]);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (unsigned int ix=0;ix<2;++ix)
        scale(_npion[ix], crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _npion[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2001_I579319);


}
