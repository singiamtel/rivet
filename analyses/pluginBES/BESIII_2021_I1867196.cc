// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> D_s* D_sJ
  class BESIII_2021_I1867196 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1867196);


    /// @name Analysis methods
    ///@{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      // Histograms
      for(unsigned int ix=0;ix<3;++ix) {
        book(_numD[ix],1+ix,1,1);
      }
      for (const string& en : _numD[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over D_s*
      for (const Particle & Dstar : ufs.particles(Cuts::abspid==433)) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(Dstar,nRes,ncount);
        bool matched=false;
        for (const Particle & p : ufs.particles(Cuts::abspid==10431 or
                                  Cuts::abspid==10433 or
                                  Cuts::abspid==20433)) {
          // check particle and antiparticle
          if(Dstar.pid()*p.pid()>0) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p,nRes2,ncount2);
          if(ncount2!=0) continue;
          matched=true;
          for (auto const & val : nRes2) {
            if(val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            if(p.abspid()==10431) _numD[0]->fill(_ecms);
            else if(p.abspid()==20433) _numD[1]->fill(_ecms);
            else if(p.abspid()==10433) _numD[2]->fill(_ecms);
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/picobarn/sumOfWeights();
      for (unsigned int ix=0;ix<3;++ix)
        scale(_numD[ix], fact);
    }

    ///@}


    /// @name Histograms
    ///@{
    BinnedHistoPtr<string> _numD[3];
    string _ecms;
    ///@}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1867196);

}
