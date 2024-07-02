// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> Ds(*) D_s1
  class BESIII_2020_I1795949 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1795949);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");

      // Book histograms
      book(_c_Ds    , 1, 1, 1);
      book(_c_DsStar, 2, 1, 1);
      for (const string& en : _c_Ds.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for(const Particle &child : p.children()) {
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
      // total analyse final state
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // unstable charm analysis
      Particles ds = apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==431 or Cuts::abspid==433 or
                                                                      Cuts::abspid==20433);
      for(unsigned int ix=0;ix<ds.size();++ix) {
        const Particle& p1 = ds[ix];
        int id1 = abs(p1.pid());
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
        bool matched=false;
        int sign = p1.pid()/id1;
        // loop over the other fs particles
        for(unsigned int iy=ix+1;iy<ds.size();++iy) {
          const Particle& p2 = ds[iy];
          fs = true;
          for (const Particle & child : p2.children()) {
            if(child.pid()==p2.pid()) {
              fs = false;
              break;
            }
          }
          if(!fs) continue;
          if(p2.pid()/abs(p2.pid())==sign) continue;
          int id2 = abs(p2.pid());
          if(!p2.parents().empty() && p2.parents()[0].pid()==p1.pid())
            continue;
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
            if((id1==431 && id2==20433) || (id1==20433 && id2==421)) {
              _c_Ds->fill(_ecms);
            }
            else if((id1==433 && id2==20433) || (id1==20433 && id2==433)) {
              _c_DsStar->fill(_ecms);
            }
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights()/picobarn;
      scale(_c_Ds    ,fact);
      scale(_c_DsStar,fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _c_Ds,_c_DsStar;
    string _ecms;
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1795949);

}
