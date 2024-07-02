// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief 2K+2K- and K+K- phi cross section
  class BESIII_2019_I1743841 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2019_I1743841);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book counters
      for(unsigned int ix=0;ix<2;++ix)
	book(_sigma[ix],ix+1, 1, 1);

      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
        if(child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else
          findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if(ntotal==4 && nCount[321]==2 && nCount[-321]==2)
        _sigma[0]->fill(_ecms);
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p :  ufs.particles(Cuts::pid==333)) {
        if(p.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // phi K+K-
        if(ncount==2) {
          bool matched = true;
          for(auto const & val : nRes) {
            if(abs(val.first)==321) {
              if(val.second!=1) {
                matched = false;
                break;
              }
            }
            else if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched)
            _sigma[1]->fill(_ecms);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
        scale(_sigma[ix],crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2019_I1743841);

}
