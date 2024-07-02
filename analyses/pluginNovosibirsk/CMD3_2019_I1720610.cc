// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> 3(pi+pi-)pi0
  class CMD3_2019_I1720610 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CMD3_2019_I1720610);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==221 || Cuts::pid==223), "UFS");

      // Book histograms
      book(_sigma[0], "TMP/c_all",   refData(1, 1, 1));
      book(_sigma[1], "TMP/c_eta",   refData(1, 1, 2));
      book(_sigma[2], "TMP/c_omega", refData(1, 1, 3));

    }

    void findChildren(const Particle& p, map<long,int>& nRes, int& ncount) {
      for (const Particle& child : p.children()) {
        if (child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else {
          findChildren(child,nRes,ncount);
        }
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      // find the final-state particles
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      if (ntotal==7 && nCount[211]==3 && nCount[-211]==3 && nCount[111] ==1 ) {
        _sigma[0]->fill(sqrtS()/MeV);
      }
      // find omega/phi + eta
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool found=false;
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // eta/omega 2(pi+pi-)
        if(ncount==4) {
          bool matched = true;
          for (const auto& val : nRes) {
            if (abs(val.first)==211 ) {
              if (val.second !=2) {
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
            if (p.pid()==221) {
              _sigma[1]->fill(sqrtS()/MeV);
            }
            else if (p.pid()==223) {
              _sigma[2]->fill(sqrtS()/MeV);
            }
            found = true;
            break;
          }
        }
        if (found) break;
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/nanobarn/sumOfWeights();
      for(unsigned int ix=0;ix<3;++ix) {
        scale(_sigma[ix], fact);
        Estimate1DPtr tmp;
        book(tmp,1,1,1+ix);
        barchart(_sigma[ix],tmp);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _sigma[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CMD3_2019_I1720610);


}
