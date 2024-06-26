// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Lambda Lambdabar cross section
  class BESIII_2018_I1627871 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1627871);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nLambda, 1, 1, 1);
      for (const string& en : _nLambda.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the Lambdas
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
        const Particle& p1 = ufs.particles()[ix];
        if(abs(p1.pid())!=3122) continue;
        bool matched = false;
        // check fs
        bool fs = true;
        for (const Particle & child : p1.children()) {
          if(child.pid()==p1.pid()) {
            fs = false;
            break;
          }
        }
        if(!fs) continue;
        // find the children
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        for(unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
          const Particle& p2 = ufs.particles()[iy];
          if(abs(p2.pid())!=3122) continue;
          // check fs
          bool fs = true;
          for (const Particle & child : p2.children()) {
            if(child.pid()==p2.pid()) {
              fs = false;
              break;
            }
          }
          if(!fs) continue;
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
            _nLambda->fill(_ecms);
            break;
          }
        }
        if(matched) break;
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


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1627871);


}
