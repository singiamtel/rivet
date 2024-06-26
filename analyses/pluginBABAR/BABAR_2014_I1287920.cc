// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+e- > K0K0 (+pions)
  class BABAR_2014_I1287920 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2014_I1287920);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<4;++ix) {
        book(_sigma[ix], 9+ix, 1, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty() && _ecms[2].empty() && _ecms[3].empty())
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

      if(ntotal==2 && nCount[130]==1 && nCount[310]==1) {
	if(!_ecms[0].empty()) _sigma[0]->fill(_ecms[0]);
      }
      else if( ntotal==4 &&  nCount[130]==1 && nCount[310]==1 && nCount[211]==1 && nCount[-211]==1) {
	if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
      }
      else if( ntotal==4 && nCount[310]==2 && nCount[211]==1 && nCount[-211]==1 ) {
	if(!_ecms[2].empty()) _sigma[2]->fill(_ecms[2]);
      }
      else if( ntotal==4 && nCount[310]==2 && nCount[321]==1 && nCount[-321]==1) {
	if(!_ecms[3].empty()) _sigma[3]->fill(_ecms[3]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<4;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[4];
    string _ecms[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2014_I1287920);


}
