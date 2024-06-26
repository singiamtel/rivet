// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > pi pi psi_2
  class BESIII_2023_I2158340 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2158340);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      //check CMS energy
      if (!inRange(sqrtS(),4.23,4.70)) {
        throw Error("Unexpected sqrtS,must be between 4.23 and 4.70 GeV");
      }
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // histos
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_n[ix], "TMP/n_"+toString(ix+1),refData(1,1,1));
      }
    }

    void findChildren(const Particle& p,map<long,int>& nRes, int& ncount) {
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
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& chi : ufs.particles(Cuts::pid==20445)) {
        if (chi.children().empty()) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(chi,nRes,ncount);
        if (ncount!=2) continue;
        if (nRes.find(211)!=nRes.end() && nRes.find(-211)!=nRes.end() && nRes[211]==1 && nRes[-211]==1) {
          _n[0]->fill(sqrtS());
          break;
        }
        else if (nRes.find(111)!=nRes.end() && nRes[111]==2) {
          _n[1]->fill(sqrtS());
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      Estimate1DPtr mult;
      book(mult, 1, 1, 1);
      divide(_n[1],_n[0],mult);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _n[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2158340);

}
