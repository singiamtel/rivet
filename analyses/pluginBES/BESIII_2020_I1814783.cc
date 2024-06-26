// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > sigma+- sigmabar -+
  class BESIII_2020_I1814783 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1814783);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_n_plus [ix],1+ix,1,1);
        book(_n_minus[ix],1+ix,1,2);
        for (const string& en : _n_plus[ix].binning().edges<0>()) {
          const double end = std::stod(en)*GeV;
          if (isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if(isCompatibleWithSqrtS(2.396*GeV)) {
        book(_h_cTheta_A,3,1,1);
        book(_h_cTheta_B,3,1,2);
      }
      if (_ecms[0].empty() && _ecms[1].empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
      }
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
      // total hadronic and muonic cross sections
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // find the Sigmas
      const FinalState& ufs = apply<UnstableParticles>(event, "UFS");
      for(unsigned int ix=0;ix<ufs.particles().size();++ix) {
        const Particle& p1 = ufs.particles()[ix];
        if(abs(p1.pid())!=3112&&abs(p1.pid())!=3222) continue;
        bool matched = false;
        // check fs
        bool fs = true;
        for(const Particle & child : p1.children()) {
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
          if(p2.pid() != -p1.pid()) continue;
          // check fs
          bool fs = true;
          for(const Particle & child : p2.children()) {
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
            if(abs(p1.pid())==3222) {
              for(unsigned int ix=0;ix<2;++ix) {
                if(!_ecms[ix].empty()) _n_plus[ix]->fill(_ecms[ix]);
              }
              if(_h_cTheta_A) {
                double cTheta = p1.pid()>0 ?
                  cos(p1.momentum().polarAngle()) :
                  cos(p2.momentum().polarAngle());
                _h_cTheta_A->fill(cTheta);
                _h_cTheta_B->fill(cTheta);
              }
            }
            else if(abs(p1.pid())==3112) {
              for(unsigned int ix=0;ix<2;++ix) {
                if(!_ecms[ix].empty()) _n_minus[ix]->fill(_ecms[ix]);
              }
            }
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      if(_h_cTheta_A) {
        normalize(_h_cTheta_A);
        normalize(_h_cTheta_B);
      }
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int iy=0;iy<2;++iy) {
        scale(_n_plus [iy],fact);
        scale(_n_minus[iy],fact);
      }
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _n_plus[2],_n_minus[2];
    Histo1DPtr _h_cTheta_A,_h_cTheta_B;
    string _ecms[2];
    /// @}

  };


  RIVET_DECLARE_PLUGIN(BESIII_2020_I1814783);

}
