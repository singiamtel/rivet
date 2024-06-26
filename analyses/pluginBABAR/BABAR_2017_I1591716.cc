// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- -> K_S0K+-pi-+ with pi0, eta 
  class BABAR_2017_I1591716 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2017_I1591716);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      
      for(unsigned int ix=0;ix<2;++ix)
        book(_sigma[ix], 1+ix, 1, 1);
      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty())
        MSG_ERROR("Beam energy incompatible with analysis.");
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
      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      // stable histos
      if (ntotal == 4 && nCount[310] == 1 && nCount[111] == 1 &&
          ( (nCount[ 321]==1 && nCount[-211]==1) ||
            (nCount[-321]==1 && nCount[ 211]==1))) {
        _sigma[0]->fill(_ecms);
      }

      // unstable particles
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if (p.children().empty()) continue;
        if (p.pid()!=221) continue;
        map<long,int> nRes=nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched  = true;
        if (p.pid()==221 && ncount==3) {
          for (const auto& val : nRes) {
            if (val.first==310 || val.first==111) {
              if (val.second!=1) {
                matched = false;
                break;
              }
            }
            else if (abs(val.first)==321 || abs(val.first)==211)
              continue;
            else if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) {
            if ((nRes[321]==1 && nRes[-211]==1 && nRes[-321]==0 && nRes[211]==0) ||
                (nRes[321]==0 && nRes[-211]==0 && nRes[-321]==1 && nRes[211]==1)) {
              _sigma[1]->fill(_ecms);
            }
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double fact = crossSection()/ sumOfWeights() /nanobarn;
      for(unsigned int ix=0;ix<2;++ix)
        scale(_sigma[ix],fact);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2017_I1591716);


}
