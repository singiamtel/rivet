// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief  e+ e- > Sigma Sigmabar
  class BELLE_2022_I2173361 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2022_I2173361);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==3222 || Cuts::abspid==3212), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_sigma[ix],"TMP/c_"+toString(ix+1),refData(ix+1, 1, 1));
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if (child.children().empty()) {
          nRes[child.pid()]-=1;
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
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
      for (unsigned int ix=0; ix<ufs.particles().size(); ++ix) {
        const Particle& p1 = ufs.particles()[ix];
        bool matched = false;
        // check fs
        bool fs = true;
        for (const Particle& child : p1.children()) {
          if (child.pid()==p1.pid()) {
            fs = false;
            break;
          }
        }
        if (!fs) continue;
        // find the children
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        for (unsigned int iy=ix+1;iy<ufs.particles().size();++iy) {
          const Particle& p2 = ufs.particles()[iy];
          if (p2.pid() != -p1.pid()) continue;
          // check fs
          bool fs = true;
          for (const Particle & child : p2.children()) {
            if (child.pid()==p2.pid()) {
              fs = false;
              break;
            }
          }
          if (!fs) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=0) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
      	  if (matched) {
      	    if(abs(p1.pid())==3212) _sigma[0]->fill(sqrtS());
            else if (abs(p1.pid())==3222) _sigma[1]->fill(sqrtS());
            break;
          }
      	}
      	if (matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int iy=0;iy<2;++iy) {
        scale(_sigma[iy], fact);
        Estimate1DPtr tmp;
        book(tmp, iy+1, 1, 1);
        barchart(_sigma[iy], tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2022_I2173361);

}
