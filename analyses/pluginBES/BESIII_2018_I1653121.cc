// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+e- > KK J/Psi
  class BESIII_2018_I1653121 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2018_I1653121);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==443), "UFS");
      book(_nKp, 1, 1, 1);
      book(_nK0, 1, 1, 2);
      for (const string& en : _nKp.binning().edges<0>()) {
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
      for (const Particle& p : fs.particles()) {
        nCount[p.pid()] += 1;
        ++ntotal;
      }
      const FinalState& ufs = apply<FinalState>(event, "UFS");
      for (const Particle& p : ufs.particles()) {
        if(p.children().empty()) continue;
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        // J/psi KK
        if(ncount!=2) continue;
        bool matched = true;
        for(auto const & val : nRes) {
          if(abs(val.first)==321 || abs(val.first)==310) {
            continue;
          }
          else if(val.second!=0) {
            matched = false;
            break;
          }
        }
        if(matched) {
          if(nRes[321]==1 && nRes[-321]==1) {
            _nKp->fill(_ecms);
            break;
          }
          else if(nRes[310]==2) {
            _nK0->fill(_ecms);
            break;
          }
        }
      }
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      const double fact = crossSection()/ sumOfWeights() /picobarn;
      scale(_nKp,fact);
      scale(_nK0,fact);
    }
    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _nKp, _nK0;
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2018_I1653121);

}
