// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

namespace Rivet {


  /// @brief e+e- -> 2pi+2pi-, K+K-pi+pi- and 2K+2K-
  class BABAR_2005_I676691 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I676691);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");

      // Book histograms
      for(unsigned int ix=0;ix<3;++ix) {
        book(_sigma[ix], 1+ix, 1, 1);
        for (const string& en : _sigma[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if (_ecms[0].empty() && _ecms[1].empty() && _ecms[2].empty())
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

      if( nCount[211]==2 && nCount[-211]==2) {
	if(!_ecms[0].empty()) _sigma[0]->fill(_ecms[0]);
      }
      else if(nCount[321]==1 && nCount[-321]==1 && nCount[211]==1 && nCount[-211]==1) {
	if(!_ecms[1].empty()) _sigma[1]->fill(_ecms[1]);
      }
      else if( nCount[321]==2 && nCount[-321]==2) {
	if(!_ecms[2].empty()) _sigma[2]->fill(_ecms[2]);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<3;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I676691);


}
