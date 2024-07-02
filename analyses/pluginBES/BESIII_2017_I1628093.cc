// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to \Lambda_c^+\bar{\Lambda}_c^-$ near threshold
  class BESIII_2017_I1628093 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2017_I1628093);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nLambda, 1, 1, 1);
      for (const string& en : _nLambda.binning().edges<0>()) {
        const double end = std::stod(en)*MeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");

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
      const FinalState& ufs = apply<FinalState>(event, "UFS");

      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
        const Particle& p1 = ufs.particles()[ix];
        if(abs(p1.pid())!=4122) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        bool matched=false;
        for(unsigned int iy=0;iy<ufs.particles().size();++iy) {
          if(ix==iy) continue;
          const Particle& p2 = ufs.particles()[iy];
          if(p2.pid() != -p1.pid()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if(ncount2!=0) continue;
          matched=true;
          for(auto const & val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            break;
          }
        }
        if(matched) {
          _nLambda->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nLambda, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nLambda;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2017_I1628093);


}
