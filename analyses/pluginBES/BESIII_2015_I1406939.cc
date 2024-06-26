// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief Cross section for $e^+e^-\to \omega\chi_{c(1,2)}$ at energies between 4.42 and 4.6 GeV
  class BESIII_2015_I1406939 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2015_I1406939);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<2;++ix) {
        book(_nChi[ix], 2+ix,1,6);
        for (const string& en : _nChi[ix].binning().edges<0>()) {
          double end = std::stod(en)*GeV;
          if(isCompatibleWithSqrtS(end)) {
            _ecms[ix] = en;
            break;
          }
        }
      }
      if(_ecms[0].empty() && _ecms[1].empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for( const Particle &child : p.children()) {
        if(child.children().empty()) {
          --nRes[child.pid()];
          --ncount;
        }
        else
          findChildren(child,nRes,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p :  fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      bool matched = false;
      for (const Particle& p :  ufs.particles(Cuts::pid==20443 or Cuts::pid==445)) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        for (const Particle& p2 :  ufs.particles(Cuts::pid==223)) {
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
            if(p.pid()==20443) {
              if(!_ecms[0].empty()) _nChi[0]->fill(_ecms[0]);
            }
            else if(p.pid()==445){
              if(!_ecms[1].empty()) _nChi[1]->fill(_ecms[1]);
            }
            break;
          }
        }
        if(matched) break;
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact =  crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<2;++ix) {
        scale(_nChi[ix],fact);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nChi[2];
    string _ecms[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2015_I1406939);


}
