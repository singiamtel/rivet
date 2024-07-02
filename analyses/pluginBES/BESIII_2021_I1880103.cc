// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableParticles.hh"

namespace Rivet {


  /// @brief e+ e- > gamma chi_c(1,2)
  class BESIII_2021_I1880103 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BESIII_2021_I1880103);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      declare(FinalState(), "FS");
      declare(UnstableParticles(Cuts::pid==20443 || Cuts::pid==445), "UFS");
      for (unsigned int ix=0; ix<2; ++ix) {
        book(_sigmaChi[ix], 1+ix, 1, 1);
      }

      for (const string& en : _sigmaChi[0].binning().edges<0>()) {
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
      for (const Particle& p : ufs.particles()) {
        if (p.children().empty()) continue;
        // find the chi_c
        map<long,int> nRes = nCount;
        int ncount = ntotal;
        findChildren(p,nRes,ncount);
        bool matched=true;
        if (ncount!=1) continue;
        for (const auto& val : nRes) {
          if (val.first==PID::PHOTON) {
            if (val.second!=1) {
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
          if (p.pid()==20443) {
            _sigmaChi[0]->fill(_ecms);
          }
          else {
            _sigmaChi[1]->fill(_ecms);
          }
          break;
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_sigmaChi, crossSection()/ sumOfWeights() /picobarn);
    }

    /// @}


    /// @name Histograms
    /// @{
    BinnedHistoPtr<string> _sigmaChi[2];
    string _ecms;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BESIII_2021_I1880103);

}
