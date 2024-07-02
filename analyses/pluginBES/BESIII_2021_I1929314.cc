// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> 4 charged particles (+pi0) cross sections
  class BESIII_2021_I1929314 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1929314);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      for (unsigned int ix=0;ix<8;++ix) {
        book(_nCharged[ix], 1, 1, ix+1);
      }

      for (const string& en : _nCharged[0].binning().edges<0>()) {
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
      if(ntotal==4) {
	if(nCount[211]==1 && nCount[-211]==1) {
	  if(nCount[321]==1 && nCount[-321]==1 )
	    _nCharged[0]->fill(_ecms);
	  else if(nCount[2212]==1 && nCount[-2212]==1 )
	    _nCharged[3]->fill(_ecms);
	}
	else if(nCount[321]==2 && nCount[-321]==2 )
	  _nCharged[1]->fill(_ecms);
	else if (nCount[211]==2 && nCount[-211]==2 )
	  _nCharged[2]->fill(_ecms);
      }
      else if(ntotal==5 && nCount[111]==1) {
	if(nCount[211]==1 && nCount[-211]==1) {
	  if(nCount[321]==1 && nCount[-321]==1 )
	    _nCharged[4]->fill(_ecms);
	  else if(nCount[2212]==1 && nCount[-2212]==1 )
	    _nCharged[7]->fill(_ecms);
	}
	else if(nCount[321]==2 && nCount[-321]==2 )
	  _nCharged[5]->fill(_ecms);
	else if (nCount[211]==2 && nCount[-211]==2 )
	  _nCharged[6]->fill(_ecms);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for (unsigned int ix=0;ix<8;++ix)
        scale(_nCharged[ix],fact);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nCharged[8];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1929314);

}
