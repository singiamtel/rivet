// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"


namespace Rivet {


  /// @brief e+ e- > K+K-, KS0 KL0 and pi+pi-pi0
  class SND_2001_I533574 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(SND_2001_I533574);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      for (unsigned int ix=0;ix<2;++ix) {
        for (unsigned int iy=0;iy<4;++iy)
          book(_sigma[ix][iy], 1+ix, 1, 1+iy);
        
        for (const string& en : _sigma[ix][3].binning().edges<0>()) {
          double end = std::stod(en)*MeV;
          if(isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
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
      if(ntotal==2) {
	if(nCount[321]==1 && nCount[-321]==1) {
          for(unsigned int ix=0;ix<2;++ix)
            _sigma[ix][0]->fill(_ecms[ix]);
        }
        else if(nCount[130]==1 && nCount[310]==1) {
          for(unsigned int ix=0;ix<2;++ix) {
            _sigma[ix][1]->fill(_ecms[ix]);
            _sigma[ix][2]->fill(_ecms[ix]);
          }
        }
      }
      else if(ntotal==3 && nCount[211] == 1 && nCount[-211] == 1 && nCount[111] == 1) {
        for(unsigned int ix=0;ix<2;++ix)
          _sigma[ix][3]->fill(_ecms[ix]);
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for (unsigned int ix=0;ix<2;++ix) {
        for (unsigned int iy=0;iy<4;++iy) {
          scale(_sigma[ix][iy],fact);
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2][4];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(SND_2001_I533574);


}
