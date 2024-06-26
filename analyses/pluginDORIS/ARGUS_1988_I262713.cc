// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> K*+K*-
  class ARGUS_1988_I262713 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(ARGUS_1988_I262713);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // book histos
      if (inRange(sqrtS()/GeV,1.4,3.0)) {
	      for (unsigned int ix=0; ix<4; ++ix)
	        book(_nMeson[ix],"TMP/nMeson_"+toString(ix+1));
      }
      else {
	      throw Error("Invalid CMS energy for ARGUS_1988_I262713");
      }
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
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
      bool resonant=false;
      // find any K* mesons
      Particles Kstar=apply<UnstableParticles>(event, "UFS").particles(Cuts::abspid==323);
      for (unsigned int ix=0;ix<Kstar.size();++ix) {
       	if(Kstar[ix].children().empty()) continue;
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(Kstar[ix],nRes,ncount);
        bool matched = false;
        for (unsigned int iy=ix+1;iy<Kstar.size();++iy) {
          if(Kstar[iy].children().empty()) continue;
          if(Kstar[ix].pid()!=-Kstar[iy].pid()) continue;
          map<long,int> nRes2=nRes;
          int ncount2 = ncount;
          findChildren(Kstar[iy],nRes2,ncount2);
          if(ncount2 !=0 ) continue;
          matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            _nMeson[0]->fill();
            resonant=true;
            break;
          }
        }
        if (matched) break;
      }
      // 4 meson final state
      if (ntotal==4) {
        if (nCount[PID::K0S]==2 &&
           nCount[PID::PIPLUS]==1 && nCount[PID::PIMINUS]==1 ) {
          _nMeson[1]->fill();
        }
        else if ( nCount[PID::K0S]==1 && nCount[PID::PI0]==1 &&
           ((nCount[PID::KPLUS ]==1 && nCount[PID::PIMINUS]==1) ||
            (nCount[PID::KMINUS]==1 && nCount[PID::PIPLUS ]==1))) {
          _nMeson[2]->fill();
          if (!resonant) _nMeson[3]->fill();
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nMeson, crossSection()/nanobarn/sumOfWeights());
      // loop over tables in paper
      for (unsigned int ix=0; ix<4; ++ix) {
        BinnedEstimatePtr<string> mult;
        book(mult, ix+1, 1, 1);
        for (auto& b : mult->bins()) {
          if (isCompatibleWithSqrtS(std::stod(b.xEdge())/GeV)) {
            b.setVal(_nMeson[ix]->val());
            b.setErr(_nMeson[ix]->err());
          }
        }
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    CounterPtr _nMeson[4];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(ARGUS_1988_I262713);

}
