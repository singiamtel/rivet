// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > eta' J/psi
  class BESIII_2020_I1762922 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2020_I1762922);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(), "UFS");
      book(_nJPsi, 1,1,1);
      for (const string& en : _nJPsi.binning().edges<0>()) {
        const double end = std::stod(en)*GeV;
        if (isCompatibleWithSqrtS(end)) {
          _ecms = en;
          break;
        }
      }
      if (_ecms.empty()) MSG_ERROR("Beam energy incompatible with analysis.");
    }

    void findChildren(const Particle & p,map<long,int> & nRes, int &ncount) {
      for (const Particle &child : p.children()) {
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
      for (const Particle& p : ufs.particles(Cuts::pid==443)) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched=false;
        for (const Particle& p2 : ufs.particles(Cuts::pid==331)) {
          if (p2.children().empty()) continue;
          bool psiParent=false;
          Particle parent=p2;
          while (!parent.parents().empty()) {
            parent=parent.parents()[0];
            if (parent.pid()==443) {
              psiParent=true;
              break;
            }
          }
          if (psiParent) continue;
          map<long,int> nRes2 = nRes;
          int ncount2 = ncount;
          findChildren(p2,nRes2,ncount2);
          if (ncount2!=0) continue;
          matched = true;
          for (const auto& val : nRes2) {
            if (val.second!=0) {
              matched = false;
              break;
            }
          }
          if(matched) {
            _nJPsi->fill(_ecms);
            break;
          }
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_nJPsi, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nJPsi;
    string _ecms;
    /// @}

  };



  RIVET_DECLARE_PLUGIN(BESIII_2020_I1762922);

}
