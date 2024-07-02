// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> rho+ rho-
  class CELLO_1989_I267081 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CELLO_1989_I267081);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if(inRange(sqrtS()/GeV,1.2,3.0)) {
        book(_nRho,"TMP/nRho");
      }
      else {
        throw Error("Invalid CMS energy for CELLO_1989_I267081");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
        if(child.children().empty()) {
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
      // find any rho mesons
      Particles rho=apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==213);
      for (unsigned int ix=0;ix<rho.size();++ix) {
       	if(rho[ix].children().empty()) continue;
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
          _nRho->fill();
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nRho, crossSection()/nanobarn/sumOfWeights());
      // loop over tables in paper
      BinnedEstimatePtr<string> mult;
      for (unsigned int ix=0; ix<2; ++ix) {
        book(mult, ix+1, 1, 1);
        for (auto& b : mult->bins()) {
          if (isCompatibleWithSqrtS(std::stod(b.xEdge()))) {
            b.set(_nRho->val(), _nRho->err());
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nRho;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CELLO_1989_I267081);

}
