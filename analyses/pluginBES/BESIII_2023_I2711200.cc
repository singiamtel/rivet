// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> pi+pi-eta
  class BESIII_2023_I2711200 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2023_I2711200);


    /// @name Analysis methods
   /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      for(unsigned int ix=0;ix<3;++ix)     
        book(_sigma[ix] , 1+ix,1,1);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        double end = std::stod(en)*GeV;
        if(isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if(_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
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
      const FinalState& fs = apply<FinalState>(event, "FS");
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      // loop over eta mesons
      for (const Particle& p : ufs.particles(Cuts::pid==221)) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched = true;
        for (const auto& val : nRes) {
          if (abs(val.first)==211) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.second!=0) {
            matched = false;
            break;
          }
        }
        if (!matched) continue;
        _sigma[0]->fill(_ecms);
        for (const Particle& p2 : ufs.particles(Cuts::pid==113)) {
          map<long,int> nResB = nRes;
          int ncountB = ncount;
          findChildren(p2,nResB,ncountB);
          if (ncountB!=0) continue;
          bool matched = true;
          for (const auto& val : nResB) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
	  if(matched) _sigma[1]->fill(_ecms);
        }
      }
      // loop over a_2 mesons
      for (const Particle& p : ufs.particles(Cuts::abspid==215 || Cuts::pid == 115)) {
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched = true;
        int ipi=111;
        if (p.abspid()==215) ipi=211*p.pid()/p.abspid();
        for (const auto& val : nRes) {
          if (abs(val.first)==ipi) {
            if (val.second !=1) {
              matched = false;
              break;
            }
          }
          else if (val.second!=0) {
            matched = false;
            break;
          }
        }
        if (matched) {
	  _sigma[2]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /picobarn;
      for(unsigned int ix=0;ix<3;++ix) {
	scale(_sigma[ix],fact);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2023_I2711200);

}
