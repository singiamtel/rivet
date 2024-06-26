// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief gamma gamma -> p K+ pbar K-
  class BELLE_2016_I1444981 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BELLE_2016_I1444981);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==102134), "UFS");
      // counters
      for (unsigned int ix=0;ix<2;++ix) {
        book(_sigma[ix],"TMP/sigma_"+toString(ix+1),refData(ix+1, 1, 1));
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
      // find the final-state particles
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // p pbar K+ K-
      if (ntotal==4 &&
          nCount[2212]==1 && nCount[-2212]==1 &&
          nCount[ 321]==1 && nCount[ -321]==1) _sigma[0]->fill(sqrtS());
      // intermediate Lambda(1520)
      for (const Particle& lam : apply<UnstableParticles>(event, "UFS").particles()) {
       	map<long,int> nRes=nCount;
       	int ncount = ntotal;
       	findChildren(lam,nRes,ncount);
        int sign = lam.pid()/lam.abspid();
        if (ncount!=2) continue;
        bool matched = true;
        for (const auto& val : nRes) {
          if (val.first==-sign*2212 || val.first==-sign*321) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
          _sigma[1]->fill(sqrtS());
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/picobarn/sumOfWeights();
      // loop over tables in paper
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_sigma[ix],fact);
        Estimate1DPtr tmp;
        book(tmp,ix+1,1,1);
        barchart(_sigma[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BELLE_2016_I1444981);

}
