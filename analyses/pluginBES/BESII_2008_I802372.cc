// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > D Dbar
  class BESII_2008_I802372 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESII_2008_I802372);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::abspid==411 || Cuts::abspid==421), "UFS");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
        book(_sigma[ix], 1, 1, 1+ix);
      }

      for (const string& en : _sigma[0].binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) {
        MSG_ERROR("Beam energy incompatible with analysis.");
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

      map<long,int> nCount;
      int ntotal(0);
      for (const Particle& p : fs.particles()) {
      	nCount[p.pid()] += 1;
      	++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle & p1 : ufs.particles()) {
        if(p1.pid()<0) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p1,nRes,ncount);
        bool matched=false;
        for (const Particle & p2 : ufs.particles()) {
          if (p2.pid()!=-p1.pid()) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=0) continue;
          matched=true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if (matched) break;
        }
        if (matched) {
          _sigma[2]->fill(_ecms);
          if (p1.abspid()==421)      _sigma[0]->fill(_ecms);
          else if (p1.abspid()==411) _sigma[1]->fill(_ecms);
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigma, crossSection()/ sumOfWeights() /nanobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigma[3];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESII_2008_I802372);

}
