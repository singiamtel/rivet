// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> rho+ rho-
  class ARGUS_1991_I315058 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1991_I315058);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,1.2,2.2)) {
        for (unsigned int ix=0; ix<3; ++ix) {
          book(_c[ix],"TMP/nMeson_"+toString(ix+1));
        }
      }
      else {
        throw Error("Invalid CMS energy for ARGUS_1991_I315058");
      }
    }

    void findChildren(const Particle& p,map<long,int> & nRes, int &ncount) {
      for (const Particle& child : p.children()) {
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      bool foundRes=false;
      // find any rho mesons
      Particles rho=apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==213);
      for (unsigned int ix=0;ix<rho.size();++ix) {
       	if (rho[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(rho[ix],nRes,ncount);
        bool matched = false;
        for (unsigned int iy=ix+1;iy<rho.size();++iy) {
          if (rho[iy].children().empty()) continue;
          if (rho[ix].pid()!=-rho[iy].pid()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(rho[iy],nRes2,ncount2);
          if (ncount2 !=0 ) continue;
          matched=true;
          for (auto const & val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            break;
          }
        }
        if (matched) {
          _c[0]->fill();
          foundRes=true;
          break;
        }
        else {
          int sign = rho[ix].pid()/rho[ix].abspid();
          bool matched2=true;
          for (const auto& val : nRes) {
            if (val.first==-sign*211 || val.first==111) {
              if (val.second!=1) {
                matched2 = false;
                break;
              }
            }
            else {
              if (val.second!=0) {
                matched2 = false;
                break;
              }
            }
          }
          if (matched2) {
            _c[1]->fill();
            foundRes=true;
            break;
          }
        }
      }
      // 4 pion final-state
      if (ntotal==4) {
        if (nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 && nCount[PID::PI0]==2) {
          if (!foundRes) _c[2]->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_c, crossSection()/nanobarn/sumOfWeights());
      // loop over tables in paper
      for (unsigned int ix=0; ix<3; ++ix) {
        Estimate1DPtr mult;
        book(mult, 1, 1, 1+ix);
        for (auto& b : mult->bins()) {
          if (inRange(sqrtS(), b.xMin(), b.xMax())) {
            b.set(_c[ix]->val(), _c[ix]->err());
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _c[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1991_I315058);

}
